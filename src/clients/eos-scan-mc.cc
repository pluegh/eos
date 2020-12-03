/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011, 2013 Danny van Dyk
 * Copyright (c) 2011 Frederik Beaujean
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <config.h>

#include <eos/constraint.hh>
#include <eos/observable.hh>
#include <eos/statistics/goodness-of-fit.hh>
#include <eos/statistics/test-statistic-impl.hh>
#include <eos/statistics/log-posterior.hh>
#include <eos/statistics/markov-chain-sampler.hh>
#if EOS_ENABLE_PMC
#include <eos/statistics/population-monte-carlo-sampler.hh>
#endif
#include <eos/utils/destringify.hh>
#include <eos/utils/instantiation_policy-impl.hh>
#include <eos/utils/hdf5.hh>
#include <eos/utils/log.hh>
#include <eos/utils/stringify.hh>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>

#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnPrint.h>

#include <time.h>

#include <gsl/gsl_cdf.h>

using namespace eos;

class DoUsage
{
    private:
        std::string _what;

    public:
        DoUsage(const std::string & what) :
            _what(what)
        {
        }

        const std::string & what() const
        {
            return _what;
        }
};

struct ObservableInput
{
        ObservablePtr observable;

        Kinematics kinematics;

        double min, central, max;
};

struct ParameterData
{
        Parameter parameter;

        double min;

        double max;

        std::string prior;
};

class CommandLine :
    public InstantiationPolicy<CommandLine, Singleton>
{
    public:
        Parameters parameters;

        Options global_options;

        LogLikelihood likelihood;

        LogPosterior log_posterior;

        MarkovChainSampler::Config mcmc_config;

#if EOS_ENABLE_PMC
        PopulationMonteCarloSampler::Config config_pmc;
#endif

        std::vector<std::shared_ptr<hdf5::File>> prerun_inputs;

        std::vector<ParameterData> scan_parameters;

        std::vector<ParameterData> nuisance_parameters;

        std::vector<ObservableInput> inputs;

        std::vector<Constraint> constraints;

        std::string creator;

        double scale_reduction;

        std::string pmc_initialization_file;
        std::string pmc_sample_file;
        bool pmc_calculate_posterior;
        unsigned pmc_calculate_posterior_min, pmc_calculate_posterior_max;
        bool pmc_draw_samples;
        bool pmc_final;
        bool pmc_update;

        bool optimize;
        std::vector<double> starting_point;

        bool goodness_of_fit;
        std::vector<double> best_fit_point;

        bool use_pmc;

        CommandLine() :
            parameters(Parameters::Defaults()),
            likelihood(parameters),
            log_posterior(likelihood),
            mcmc_config(MarkovChainSampler::Config::Quick()),
#if EOS_ENABLE_PMC
            config_pmc(PopulationMonteCarloSampler::Config::Default()),
#endif
            scale_reduction(1),
            pmc_calculate_posterior(false),
            pmc_calculate_posterior_min(0),
            pmc_calculate_posterior_max(0),
            pmc_draw_samples(false),
            pmc_final(false),
            pmc_update(false),
            optimize(false),
            goodness_of_fit(false),
            use_pmc(false)
        {
            mcmc_config.number_of_chains = 4;
            mcmc_config.need_prerun = true;
            mcmc_config.chunk_size = 1000;
            mcmc_config.parallelize = false;
            mcmc_config.use_strict_rvalue_definition = true;
        }

        void parse(int argc, char ** argv)
        {
            Log::instance()->set_log_level(ll_informational);
            Log::instance()->set_program_name("eos-scan-mc");

            std::shared_ptr<Kinematics> kinematics(new Kinematics);

            creator = std::string(argv[0]);
            for (int i = 1 ; i < argc ; ++i)
            {
                creator += ' ' + std::string(argv[i]);
            }

            for (char ** a(argv + 1), **a_end(argv + argc); a != a_end; ++a)
            {
                std::string argument(*a);

                /*
                 * format: N_SIGMAS in [0, 10]
                 * a) --scan PAR N_SIGMAS --prior ...
                 * b) --scan PAR MIN MAX  --prior ...
                 * c) --scan PAR HARD_MIN HARD_MAX N_SIGMAS --prior ...
                 */
                if (("--scan" == argument) || ("--nuisance" == argument))
                {
                    std::string name = std::string(*(++a));

                    double min = -std::numeric_limits<double>::max();
                    double max =  std::numeric_limits<double>::max();

                    // first word has to be a number
                    double number = destringify<double>(*(++a));

                    std::string keyword = std::string(*(++a));

                    VerifiedRange<double> n_sigmas(0, 10, 0);

                    // case a)
                    if ("--prior" == keyword)
                    {
                        n_sigmas = VerifiedRange<double>(0, 10, number);
                        if (n_sigmas == 0)
                            throw DoUsage("number of sigmas: number expected");
                    }
                    else
                    {
                        // case b), c)
                        min = number;
                        max = destringify<double>(keyword);

                        keyword = std::string(*(++a));

                        // watch for case c)
                        if ("--prior" != keyword)
                        {
                            n_sigmas = VerifiedRange<double>(0, 10,  destringify<double>(keyword));
                            if (n_sigmas == 0)
                                throw DoUsage("number of sigmas: number expected");
                            keyword = std::string(*(++a));
                        }
                    }

                    if ("--prior" != keyword)
                        throw DoUsage("Missing correct prior specification for '" + name + "'!");

                    std::string prior_type = std::string(*(++a));

                    LogPriorPtr prior;

                    ParameterRange range{ min, max };

                    if (prior_type == "gaussian")
                    {
                        double lower = destringify<double> (*(++a));
                        double central = destringify<double> (*(++a));
                        double upper = destringify<double> (*(++a));

                        // adjust range, but always stay within hard bound supplied by the user
                        if (n_sigmas > 0)
                        {
                            range.min = std::max(range.min, central - n_sigmas * (central - lower));
                            range.max = std::min(range.max, central + n_sigmas * (upper - central));
                        }
                        if (prior_type == "gaussian")
                        {
                            prior = LogPrior::Gauss(parameters, name, range, lower, central, upper);
                        }
                    }
                    else if (prior_type == "flat")
                    {
                        if (n_sigmas > 0)
                            throw DoUsage("Can't specify number of sigmas for flat prior");
                        prior = LogPrior::Flat(parameters, name, range);
                    }
                    else
                    {
                        throw DoUsage("Unknown prior distribution: " + prior_type);
                    }

                    bool nuisance = ("--nuisance" == argument) ? true : false;

                    if (nuisance)
                    {
                        nuisance_parameters.push_back(ParameterData{ parameters[name], range.min, range.max, prior_type });
                    }
                    else
                    {
                        scan_parameters.push_back(ParameterData{ parameters[name], range.min, range.max, prior_type });
                    }

                    // check for error in setting the prior and adding the parameter
                    if (! log_posterior.add(prior, nuisance))
                        throw DoUsage("Error in assigning " + prior_type + " prior distribution to '" + name +
                                      "'. Perhaps '" + name + "' appears twice in the list of parameters?");

                    continue;
                }

                if ("--chains" == argument)
                {
                    mcmc_config.number_of_chains = destringify<unsigned>(*(++a));
                    continue;
                }

                if ("--chunk-size" == argument)
                {
                    mcmc_config.chunk_size = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--chunks" == argument)
                {
                    mcmc_config.chunks = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--constraint" == argument)
                {
                    std::string constraint_name(*(++a));

                    Constraint c(Constraint::make(constraint_name, global_options));
                    likelihood.add(c);
                    constraints.push_back(c);

                    continue;
                }

                if ("--debug" == argument)
                {
                    Log::instance()->set_log_level(ll_debug);

                    continue;
                }

                if ("--fix" == argument)
                {
                    std::string par_name = std::string(*(++a));
                    double value = destringify<double> (*(++a));
                    log_posterior.parameters()[par_name]=value;

                    continue;
                }

                if ("--kinematics" == argument)
                {
                    std::string name = std::string(*(++a));
                    double value = destringify<double> (*(++a));
                    kinematics->declare(name);
                    kinematics->set(name, value);

                    continue;
                }

                if ("--global-option" == argument)
                {
                    std::string name(*(++a));
                    std::string value(*(++a));

                    if (! constraints.empty())
                    {
                        Log::instance()->message("eos-scan-mc", ll_warning)
                            << "Global option (" << name << " = " << value <<") only applies to observables/constraints defined from now on, "
                            << "but doesn't affect the " << constraints.size() << " previously defined constraints.";
                    }

                    global_options.set(name, value);

                    continue;
                }

                if ("--goodness-of-fit" == argument)
                {
                    // best-fit point is optional
                    goodness_of_fit = true;

                    std::string lbrace(*(++a));
                    if ("{" != lbrace)
                    {
                        --a;
                        continue;
                    }

                    // parse starting point
                    do
                    {
                        std::string word(*(++a));
                        if ("}" == word)
                            break;

                        double value = destringify<double>(word);
                        best_fit_point.push_back(value);
                    }
                    while (true);

                    continue;
                }

#if EOS_ENABLE_PMC
                if ("--hc-patch-length" == argument)
                {
                	config_pmc.patch_length = destringify<double> (*(++a));

                	continue;
                }

                if ("--hc-skip-initial" == argument)
                {
                	config_pmc.skip_initial = destringify<double> (*(++a));

                	continue;
                }

                if ("--hc-target-ncomponents" == argument)
                {
                    config_pmc.target_ncomponents = destringify<unsigned>(*(++a));

                    continue;
                }
#endif

                if ("--no-prerun" == argument)
                {
                    mcmc_config.need_prerun = false;

                    continue;
                }

                if ("--observable" == argument)
                {
                    std::string observable_name(*(++a));

                    ObservableInput input;
                    input.kinematics = *kinematics;
                    input.observable = Observable::make(observable_name, parameters,
                            *kinematics, global_options);
                    if (!input.observable)
                        throw DoUsage("Unknown observable '" + observable_name + "'");

                    input.min = destringify<double> (*(++a));
                    input.central = destringify<double> (*(++a));
                    input.max = destringify<double> (*(++a));

                    likelihood.add(input.observable, input.min, input.central, input.max);

                    inputs.push_back(input);
                    kinematics.reset(new Kinematics);

                    continue;
                }

                if ("--observable-prior" == argument)
                {
                    std::string observable_name(*(++a));

                    ObservableInput input;
                    input.kinematics = *kinematics;
                    input.observable = Observable::make(observable_name, parameters,
                            *kinematics, global_options);
                    if (!input.observable)
                        throw DoUsage("Unknown observable '" + observable_name + "'");

                    input.min = destringify<double> (*(++a));
                    input.central = destringify<double> (*(++a));
                    input.max = destringify<double> (*(++a));

                    // TODO: Hack! This is only used for putting parts of the
                    // prior into the likelihood for correlated prior information.
                    likelihood.add(input.observable, input.min, input.central, input.max, 0);

                    inputs.push_back(input);
                    kinematics.reset(new Kinematics);

                    continue;
                }

                if ("--optimize" == argument)
                {
                    optimize = true;

                    // starting point is optional
                    std::string lbrace(*(++a));
                    if ("{" != lbrace)
                    {
                        --a;
                        continue;
                    }

                    // parse starting point
                    do
                    {
                        std::string word(*(++a));
                        if ("}" == word)
                            break;

                        double value = destringify<double>(word);
                        starting_point.push_back(value);
                    }
                    while (true);

                    continue;
                }

                if ("--output" == argument)
                {
                    std::string filename(*(++a));
                    mcmc_config.output_file = filename;
#if EOS_ENABLE_PMC
                    config_pmc.output_file = filename;
#endif

                    continue;
                }

                if ("--parallel" == argument)
                {
                    mcmc_config.parallelize = false;
#if EOS_ENABLE_PMC
                    config_pmc.parallelize = mcmc_config.parallelize;
#endif

                    continue;
                }

                if ("--use-pmc" == argument)
                {
                    use_pmc = true;

                    continue;
                }

#if EOS_ENABLE_PMC
                if ("--pmc-adjust-sample-size" == argument)
                {
                    config_pmc.adjust_sample_size = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--pmc-crop-highest-weights" == argument)
                {
                    config_pmc.crop_highest_weights = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--pmc-dof" == argument)
                {
                    config_pmc.degrees_of_freedom = destringify<long>(*(++a));

                    continue;
                }

                if ("--pmc-calculate-posterior" == argument)
                {
                    pmc_calculate_posterior = true;

                    // read samples from this file
                    pmc_sample_file = std::string(*(++a));
                    pmc_calculate_posterior_min = destringify<unsigned>(*(++a));
                    pmc_calculate_posterior_max = destringify<unsigned>(*(++a));

                    // read components from this file
                    pmc_initialization_file = pmc_sample_file;

                    continue;
                }

                if ("--pmc-draw-samples" == argument)
                {
                    // samples are to be stored in ordinary output file via config
                    pmc_draw_samples = true;

                    continue;
                }

                if ("--pmc-final" == argument)
                {
                    pmc_final = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--pmc-ignore-group" == argument)
                {
                    config_pmc.ignore_groups.push_back(destringify<unsigned>(*(++a)));

                    continue;
                }

                if ("--pmc-initialize-from-file" == argument)
                {
                    pmc_initialization_file = std::string(*(++a));

                    continue;
                }

                if ("--pmc-group-by-r-value" == argument)
                {
                    config_pmc.group_by_r_value = destringify<double>(*(++a));

                    continue;
                }

                if ("--pmc-r-value-no-nuisance" == argument)
                {
                    config_pmc.r_value_no_nuisance = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--pmc-final-samples" == argument)
                {
                    config_pmc.final_samples = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--pmc-ignore-ess" == argument)
                {
                    config_pmc.ignore_eff_sample_size = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--pmc-max-updates" == argument)
                {
                    config_pmc.max_updates = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--pmc-relative-std-deviation-over-last-steps" == argument)
                {
                    config_pmc.maximum_relative_std_deviation = destringify<double>(*(++a));
                    config_pmc.minimum_steps = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--pmc-samples-per-component" == argument)
                {
                    config_pmc.samples_per_component = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--pmc-update" == argument)
                {
                    pmc_update = true;
                    pmc_initialization_file = std::string(*(++a));

                    continue;
                }
#endif
                // todo rename here and in scripts
                if ("--prerun-chains-per-partition" == argument)
                {
                    mcmc_config.number_of_chains = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--prerun-max" == argument)
                {
                    mcmc_config.prerun_iterations_max = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--prerun-min" == argument)
                {
                    mcmc_config.prerun_iterations_min = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--prerun-only" == argument)
                {
                    mcmc_config.need_prerun = true;
                    mcmc_config.store_prerun = true;
                    mcmc_config.need_main_run = false;

                    continue;
                }

                if ("--prerun-update" == argument)
                {
                    mcmc_config.prerun_iterations_update = destringify<unsigned>(*(++a));

                    continue;
                }

                if ("--print-args" == argument)
                {
                    // print arguments and quit
                   for (int i = 1 ; i < argc ; i++)
                    {
                       std::cout << "'" << argv[i] << "' ";
                    }

                    std::cout << std::endl;
                    abort();

                    continue;
                }

                if ("--proposal" == argument)
                {
                    mcmc_config.proposal = *(++a);

                    if (mcmc_config.proposal == "MultivariateStudentT")
                    {
                        double dof = destringify<double>(*(++a));
                        if (dof <= 0)
                        {
                            throw DoUsage("No (or non-positive) degree of freedom for MultivariateStudentT specified");
                        }
                        mcmc_config.student_t_degrees_of_freedom = dof;
                    }

                    continue;
                }

                if ("--seed" == argument)
                {
                    std::string value(*(++a));

                    if ("time" == value)
                    {
                        mcmc_config.seed = ::time(0);
#if EOS_ENABLE_PMC
                        config_pmc.seed = ::time(0);
#endif
                    }
                    else
                    {
                        mcmc_config.seed = destringify<unsigned long>(value);
#if EOS_ENABLE_PMC
                        config_pmc.seed = destringify<unsigned long>(value);
#endif
                    }

                    continue;
                }

                if ("--scale-reduction" == argument)
                {
                    scale_reduction = destringify<double>(*(++a));

                    continue;
                }

                if ("--store-prerun" == argument)
                {
                    mcmc_config.store_prerun = true;

                    continue;
                }

                throw DoUsage("Unknown command line argument: " + argument);
            }
        }
};

int main(int argc, char * argv[])
{
    try
    {
        auto inst = CommandLine::instance();
        inst->parse(argc, argv);

        if (inst->inputs.empty() && inst->constraints.empty())
            throw DoUsage("Neither inputs nor constraints specified");

        if (inst->nuisance_parameters.empty() &&
            inst->scan_parameters.empty())
           throw  DoUsage("Neither scan nor nuisance parameters defined");

        std::cout << std::scientific;
        std::cout << "# Scan generated by eos-scan-mc" << std::endl;
        if ( ! inst->scan_parameters.empty())
        {
            std::cout << "# Scan parameters (" << inst->scan_parameters.size() << "):" << std::endl;
            for (auto d = inst->log_posterior.parameter_descriptions().cbegin(), d_end = inst->log_posterior.parameter_descriptions().cend() ;
                 d != d_end ; ++d)
            {
                if (d->nuisance)
                    continue;
                std::cout << "#   " << inst->log_posterior.log_prior(d->parameter->name())->as_string() << std::endl;
            }
        }

        if ( ! inst->nuisance_parameters.empty())
        {
            std::cout << "# Nuisance parameters (" << inst->nuisance_parameters.size() << "):" << std::endl;
            for (auto d = inst->log_posterior.parameter_descriptions().cbegin(), d_end = inst->log_posterior.parameter_descriptions().cend() ;
                 d != d_end ; ++d)
            {
                if ( ! d->nuisance)
                    continue;
                std::cout << "#   " << inst->log_posterior.log_prior(d->parameter->name())->as_string() << std::endl;
            }
        }

        if ( ! inst->inputs.empty())
        {
            std::cout << "# Manual inputs (" << inst->inputs.size() << "):" << std::endl;
            for (auto i = inst->inputs.cbegin(), i_end = inst->inputs.cend() ; i != i_end ; ++i)
            {
                std::cout << "#   " << i->observable->name() << '['
                    << i->kinematics.as_string() << "] = (" << i->min << ", "
                    << i->central << ", " << i->max << ')' << std::endl;
            }
        }

        if ( ! inst->constraints.empty())
        {
            std::cout << "# Constraints (" << inst->constraints.size() << "):" << std::endl;
            for (auto c = inst->constraints.cbegin(), c_end = inst->constraints.cend() ; c != c_end ; ++c)
            {
                std::cout << "#  " << c->name() << ": ";
                for (auto o = c->begin_observables(), o_end = c->end_observables(); o != o_end ; ++o)
                {
                    std::cout << (**o).name() << '['
                        << (**o).kinematics().as_string() << ']'
                        << " with options: " << (**o).options().as_string();
                }
                for (auto b = c->begin_blocks(), b_end = c->end_blocks(); b != b_end ; ++b)
                {
                    std::cout << ", " << (**b).as_string();
                }
                std::cout << std::endl;
            }
        }

        // run optimization. Use starting point if given, else sample a point from the prior.
        // Optionally calculate a p-value at the mode.
        if (inst->optimize)
        {
            LogPosterior ana = inst->log_posterior;
            if (inst->starting_point.empty())
            {
                gsl_rng * rng = gsl_rng_alloc(gsl_rng_mt19937);
                gsl_rng_set(rng, ::time(0));
                for (auto i = ana.parameter_descriptions().begin(), i_end = ana.parameter_descriptions().end() ; i != i_end ; ++i)
                {
                    LogPriorPtr prior = ana.log_prior(i->parameter->name());
                    inst->starting_point.push_back(prior->sample(rng));
                }
            }

            if (inst->starting_point.size() != ana.parameter_descriptions().size())
            {
                throw DoUsage("Starting point size of" + stringify(inst->starting_point.size())
                              + " doesn't match with analysis size of " + stringify(ana.parameter_descriptions().size()));
            }

            std::cout << std::endl;
            std::cout << "# Starting optimization at " << stringify_container(inst->starting_point, 4) << std::endl;
            std::cout << std::endl;

            auto options = LogPosterior::OptimizationOptions::Defaults();
            auto ret = ana.optimize_minuit(inst->starting_point, options);

            Log::instance()->message("eos-scan-mc", ll_informational)
                << "Result from minuit:" << ret << ret.UserCovariance();
            if (inst->goodness_of_fit && inst->best_fit_point.empty())
                ana.goodness_of_fit(ret.UserParameters().Params(), 1e5, inst->mcmc_config.output_file);

            Log::instance()->message("eos-scan-mc", ll_informational)
                << "Best result: log(posterior) at "
                << stringify_container(ret.UserParameters().Params(), 6)
                << " = " << -1.0 * ret.Fval();

            {
                auto gof = GoodnessOfFit(CommandLine::instance()->log_posterior);
                Log::instance()->message("eos-scan-mc", ll_informational)
                    << "Goodness-of-Fit summary";
                for (auto c = gof.begin_chi_square(), c_end = gof.end_chi_square() ; c != c_end ; ++c)
                {
                    Log::instance()->message("eos-scan-mc", ll_informational)
                        << "  " << c->first << " : chi^2 = " << c->second.chi2;
                }
                Log::instance()->message("eos-scan-mc", ll_informational)
                    << "----------------";
                Log::instance()->message("eos-scan-mc", ll_informational)
                    << " total chi^2  = " << gof.total_chi_square();
                Log::instance()->message("eos-scan-mc", ll_informational)
                    << " total d.o.f. = " << gof.total_degrees_of_freedom();
                Log::instance()->message("eos-scan-mc", ll_informational)
                    << " p value      = " << gsl_cdf_chisq_Q(gof.total_chi_square(), gof.total_degrees_of_freedom());
            }

            return EXIT_SUCCESS;
        }

        // goodness-of-fit for user specified parameter point
        if (inst->goodness_of_fit)
        {
            inst->log_posterior.goodness_of_fit(inst->best_fit_point, 1e5, inst->mcmc_config.output_file);
            auto gof = GoodnessOfFit(CommandLine::instance()->log_posterior);
            Log::instance()->message("eos-scan-mc", ll_informational)
                << "Goodness-of-Fit summary";
            for (auto c = gof.begin_chi_square(), c_end = gof.end_chi_square() ; c != c_end ; ++c)
            {
                Log::instance()->message("eos-scan-mc", ll_informational)
                    << "  " << c->first << " : chi^2 = " << c->second.chi2;
            }
            Log::instance()->message("eos-scan-mc", ll_informational)
                << "----------------";
            Log::instance()->message("eos-scan-mc", ll_informational)
                << " total chi^2  = " << gof.total_chi_square();
            Log::instance()->message("eos-scan-mc", ll_informational)
                << " total d.o.f. = " << gof.total_degrees_of_freedom();
            Log::instance()->message("eos-scan-mc", ll_informational)
                << " p value      = " << gsl_cdf_chisq_Q(gof.total_chi_square(), gof.total_degrees_of_freedom());

            return EXIT_SUCCESS;
        }

#if EOS_ENABLE_PMC
        if (CommandLine::instance()->use_pmc)
        {
            PopulationMonteCarloSampler pop_sampler(inst->log_posterior.clone(), hdf5::File::Open(inst->pmc_initialization_file),
            		inst->config_pmc, inst->pmc_update);

            if (inst->pmc_final)
            {
                PopulationMonteCarloSampler::Status status =  pop_sampler.status();
                status.converged = true;
                pop_sampler.status(status);
            }
            if (inst->pmc_draw_samples)
            {
                pop_sampler.draw_samples();
            }
            else if (inst->pmc_calculate_posterior)
            {
                pop_sampler.calculate_weights(inst->pmc_sample_file,
                    inst->pmc_calculate_posterior_min, inst->pmc_calculate_posterior_max);
            }
            else if (inst->pmc_update)
                return EXIT_SUCCESS;
            else
                pop_sampler.run();

            return EXIT_SUCCESS;
        }
#endif

        MarkovChainSampler sampler(inst->log_posterior.clone(), inst->mcmc_config);

        sampler.run();
    }
    catch (DoUsage & e)
    {
    	// todo update
        std::cout << e.what() << std::endl;
        std::cout << "Usage: eos-scan-mc" << std::endl;
        std::cout << "  [ [--kinematics NAME VALUE]* --observable NAME LOWER CENTRAL UPPER]+" << std::endl;
        std::cout << "  [--constraint NAME]+" << std::endl;
        std::cout << "  [ [ [--scan PARAMETER MIN MAX] | [--nuisance PARAMETER MIN MAX] ] --prior [flat | [gaussian LOWER CENTRAL UPPER] ] ]+" << std::endl;
        std::cout << "  [--chains VALUE]" << std::endl;
        std::cout << "  [--chunks VALUE]" << std::endl;
        std::cout << "  [--chunksize VALUE]" << std::endl;
        std::cout << "  [--debug]" << std::endl;
        std::cout << "  [--fix PARAMETER VALUE]+" << std::endl;
        std::cout << "  [--goodness_of_fit [{ PAR_VALUE1 PAR_VALUE2 ... PAR_VALUEN }]]" << std::endl;
        std::cout << "  [--no-prerun]" << std::endl;
        std::cout << "  [--optimize [{ PAR_VALUE1 PAR_VALUE2 ... PAR_VALUEN }]]" << std::endl;
        std::cout << "  [--output FILENAME]" << std::endl;
        std::cout << "  [--scale VALUE]" << std::endl;
        std::cout << "  [--seed LONG_VALUE]" << std::endl;
        std::cout << "  [--store-prerun]" << std::endl;


        std::cout << std::endl;
        std::cout << "Example:" << std::endl;
        std::cout << "  eos-scan-mc --kinematics s_min 14.18 --kinematics s_max 16.00 \\" << std::endl;
        std::cout << "      --observable \"B->K^*ll::BR@LowRecoil\" 0.5e-7 1.25e-7 2.0e-7 \\" << std::endl;
        std::cout << "      --constraint \"B^0->K^*0gamma::BR@BaBar-2009\" \\" << std::endl;
        std::cout << "      --scan     \"Abs{c9}\"        0.0 15.0     --prior flat\\" << std::endl;
        std::cout << "      --scan     \"Arg{c9}\"        0.0  6.28319 --prior flat\\" << std::endl;
        std::cout << "      --nuisance \"mass::b(MSbar)\" 3.8  5.0     --prior gaussian 4.14 4.27 4.37" << std::endl;

        return EXIT_FAILURE;
    }
    catch (Exception & e)
    {
        std::cerr << "Caught exception: '" << e.what() << "'" << std::endl;
        return EXIT_FAILURE;
    }
    catch (...)
    {
        std::cerr << "Aborting after unknown exception" << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
