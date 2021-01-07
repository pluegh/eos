/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2018 Danny van Dyk
 * Copyright (c) 2018 Nico Gubernari
 * Copyright (c) 2018 Ahmet Kokulu
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_ANALYTIC_B_TO_P_LCSR_HH
#define EOS_GUARD_EOS_FORM_FACTORS_ANALYTIC_B_TO_P_LCSR_HH 1

#include <eos/form-factors/mesonic.hh>
#include <eos/utils/diagnostics.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/options.hh>

namespace eos
{
    namespace lcsr
    {
        struct BToPi;
        struct BToK;
        struct BToD;
        struct BsToK;
        struct BsToDs;
        struct BcToD;
    }

    template <typename Process_>
    class AnalyticFormFactorBToPLCSR :
        public FormFactors<PToP>,
        PrivateImplementationPattern<AnalyticFormFactorBToPLCSR<Process_>>
    {
        public:
            AnalyticFormFactorBToPLCSR(const Parameters &, const Options &);

            ~AnalyticFormFactorBToPLCSR();

            static FormFactors<PToP> * make(const Parameters &, const Options &);

            /* Form factors */
            virtual double f_p(const double & q2) const;
            virtual double f_0(const double & q2) const;
            virtual double f_t(const double & q2) const;
            virtual double f_m(const double & q2) const;

            /* First moments of the sum rules */
            double normalized_moment_1_f_p(const double & q2) const;
            double normalized_moment_1_f_pm(const double & q2) const;
            double normalized_moment_1_f_t(const double & q2) const;

            /* Diagnostics for unit tests */
            Diagnostics diagnostics() const;
    };

    extern template class AnalyticFormFactorBToPLCSR<lcsr::BToPi>;
    extern template class AnalyticFormFactorBToPLCSR<lcsr::BToK>;
    extern template class AnalyticFormFactorBToPLCSR<lcsr::BToD>;
    extern template class AnalyticFormFactorBToPLCSR<lcsr::BsToK>;
    extern template class AnalyticFormFactorBToPLCSR<lcsr::BsToDs>;
    extern template class AnalyticFormFactorBToPLCSR<lcsr::BcToD>;
}
#endif
