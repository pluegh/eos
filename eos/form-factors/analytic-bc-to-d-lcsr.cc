/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2020 Danny van Dyk
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

#include <eos/form-factors/analytic-b-to-p-lcsr-impl.hh>

namespace eos
{
    namespace lcsr
    {
        struct BcToD
        {
            constexpr static const char * B    = "B_c";
            constexpr static const char * m_B  = "mass::B_c";
            constexpr static const char * f_B  = "decay-constant::B_c";
            constexpr static const char * P    = "D";
            constexpr static const char * m_P  = "mass::D^0";
            constexpr static const char * f_P  = "decay-constant::D_u";
            constexpr static const char   q_v  = 'u';
            constexpr static const char   q_s  = 'c';
            constexpr static const double chi2 = 1.0;
        };

        // Bc -> D
        constexpr const char * BcToD::B;
        constexpr const char * BcToD::m_B;
        constexpr const char * BcToD::f_B;
        constexpr const char * BcToD::P;
        constexpr const char * BcToD::m_P;
        constexpr const char * BcToD::f_P;
        constexpr const char   BcToD::q_s;
    }

    template class AnalyticFormFactorBToPLCSR<lcsr::BcToD>;
}
