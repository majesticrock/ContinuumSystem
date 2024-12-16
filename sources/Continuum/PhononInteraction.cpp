#include "PhononInteraction.hpp"
#include "SCModel.hpp"
#include <Utility/Numerics/Roots/Bisection.hpp>
#include <Utility/Numerics/Interpolation.hpp>
#include <cmath>

namespace Continuum {
    void PhononInteraction::compute_renormalization_table()
    {
        assert(parent != nullptr);
        renormalization_cache.resize(parent->momentumRanges.size(), std::array<c_float, 3>{ c_float{}, c_float{}, c_float{} });

#ifdef PHONON_SC_CHANNEL_ONLY
        return;
#else
        for (MomentumIterator it(& parent->momentumRanges); it < MomentumIterator::max_idx(); ++it) {
            auto integrand_infinity = [this, &it](c_float q) -> c_float {
		    	return -q * q / std::copysign(std::abs(alpha_CUT(q, it.k)) + CUT_REGULARIZATION, alpha_CUT(q, it.k)) ;
		    };

            auto integrand_fock = [this, &it](c_float q) -> c_float {
		    	return q * q * 
		    		( 1. / std::copysign(std::abs(beta_CUT(q, it.k)) + CUT_REGULARIZATION, beta_CUT(q, it.k)) 
		    		- 1. / std::copysign(std::abs(alpha_CUT(q, it.k)) + CUT_REGULARIZATION, alpha_CUT(q, it.k)) );
		    };

            const auto singularities = get_singularities(it.k);
            // Fock part
            if (singularities[0] < parent->fermi_wavevector) {
                renormalization_cache[it.idx][0] = parent->momentumRanges.integrate(integrand_fock, parent->momentumRanges.K_MIN, singularities[0]);

                if (singularities[1] < parent->fermi_wavevector) {
                    renormalization_cache[it.idx][0] += parent->momentumRanges.integrate(integrand_fock, singularities[0], singularities[1]) +
                            parent->momentumRanges.integrate(integrand_fock, singularities[1], parent->fermi_wavevector);
                }
                else {
                    renormalization_cache[it.idx][0] += parent->momentumRanges.integrate(integrand_fock, singularities[0], parent->fermi_wavevector);
                }
            }
            else {
                renormalization_cache[it.idx][0] = parent->momentumRanges.integrate(integrand_fock, parent->momentumRanges.K_MIN, parent->fermi_wavevector);
            }

            // CUT flow part
            renormalization_cache[it.idx][1] = parent->momentumRanges.integrate(integrand_infinity, parent->momentumRanges.K_MIN, singularities[0])
                + parent->momentumRanges.integrate(integrand_infinity, singularities[0], singularities[1])
                + parent->momentumRanges.integrate(integrand_infinity, singularities[1], parent->momentumRanges.K_MAX);

        }
#endif
    }

    void PhononInteraction::compute_singularities()
    {
        assert(parent != nullptr);
        singularities_cache.resize(parent->momentumRanges.size(), std::array<c_float, 2>{ c_float{}, c_float{} });
        for (MomentumIterator it(& parent->momentumRanges); it < MomentumIterator::max_idx(); ++it) 
        {
            try {
			    singularities_cache[it.idx][0] = Utility::Numerics::Roots::bisection([this, &it](c_float q) { return beta_CUT(q, it.k);  }, parent->momentumRanges.K_MIN, parent->momentumRanges.K_MAX, PRECISION, 200);
            } catch (Utility::Numerics::Roots::NoRootException const & e) {
                singularities_cache[it.idx][0] = 2 * parent->momentumRanges.K_MAX;
            }
            try {
                singularities_cache[it.idx][1] = Utility::Numerics::Roots::bisection([this, &it](c_float q) { return alpha_CUT(q, it.k); }, parent->momentumRanges.K_MIN, parent->momentumRanges.K_MAX, PRECISION, 200);
            } catch (Utility::Numerics::Roots::NoRootException const & e) {
                singularities_cache[it.idx][1] = 2 * parent->momentumRanges.K_MAX;
            }

            // Assure that the singularities are ordered
            if (singularities_cache[it.idx][0] > singularities_cache[it.idx][1]) {
                std::swap(singularities_cache[it.idx][0], singularities_cache[it.idx][1]);
            }
		};
    }

    std::array<c_float, 2> PhononInteraction::get_singularities(c_float k) const
    {
        assert(parent != nullptr);
        int index = parent->momentumRanges.momentum_to_floor_index(k);
        if(index >= parent->momentumRanges.size() - 1) {
            index = parent->momentumRanges.size() - 2;
        }
        return {
            Utility::Numerics::linearly_interpolate(k, parent->momentumRanges[index], parent->momentumRanges[index + 1], 
                                singularities_cache[index][0], singularities_cache[index + 1][0]),
            Utility::Numerics::linearly_interpolate(k, parent->momentumRanges[index], parent->momentumRanges[index + 1], 
                                singularities_cache[index][1], singularities_cache[index + 1][1])
        };
    }

    void PhononInteraction::set_parent(SCModel const *_parent)
    {
        parent = _parent;
        ptr_rho_F = &_parent->rho_F;
        ptr_fermi_wavevector = &_parent->fermi_wavevector;
        ptr_omega_debye = &_parent->omega_debye;
        ptr_phonon_coupling = &_parent->phonon_coupling;
        ptr_momentumRanges = &_parent->momentumRanges;

        this->compute_singularities();
        this->compute_renormalization_table();
    }

    c_float PhononInteraction::alpha_CUT(c_float k, c_float k_prime) const {
        assert(parent != nullptr);
		return parent->delta_epsilon(k, k_prime) + parent->omega_debye;
	}
	c_float PhononInteraction::beta_CUT(c_float k, c_float k_prime) const {
        assert(parent != nullptr);
		return parent->delta_epsilon(k, k_prime) - parent->omega_debye;
	}

    c_float PhononInteraction::sc_channel_lower_bound(c_float k) const
	{
        assert(parent != nullptr);
#ifdef approximate_theta
		const c_float ALPHA = 2. * parent->fermi_energy - 2. * omega_debye;
#else
		const c_float ALPHA = parent->phonon_boundary_a(k) - 2. * parent->omega_debye;
#endif
		auto func = [&](c_float l) {
			return parent->phonon_boundary_b(l, ALPHA);
			};
#ifdef approximate_theta
		return Utility::Numerics::Roots::bisection(func, parent->momentumRanges.K_MIN, fermi_wavevector, PRECISION, 200);
#else
		const auto lb = func(parent->momentumRanges.K_MIN);
		const auto ub = func(k);
		if (lb * ub >= c_float{}) return parent->momentumRanges.K_MIN;
		return Utility::Numerics::Roots::bisection(func, parent->momentumRanges.K_MIN, k, PRECISION, 200);
#endif
	}

	c_float PhononInteraction::sc_channel_upper_bound(c_float k) const
	{
        assert(parent != nullptr);
#ifdef approximate_theta
		const c_float ALPHA = 2. * parent->fermi_energy + 2. * omega_debye;
#else
		const c_float ALPHA = parent->phonon_boundary_a(k) + 2. * parent->omega_debye;
#endif
		auto func = [&](c_float l) {
			return parent->phonon_boundary_b(l, ALPHA);
			};
#ifdef approximate_theta
		return Utility::Numerics::Roots::bisection(func, fermi_wavevector, parent->momentumRanges.K_MAX, PRECISION, 200);
#else
		const auto lb = func(k);
		const auto ub = func(parent->momentumRanges.K_MAX);
		if (lb * ub >= c_float{}) return parent->momentumRanges.K_MAX;
		return Utility::Numerics::Roots::bisection(func, k, parent->momentumRanges.K_MAX, PRECISION, 200);
#endif
	}

    c_float PhononInteraction::renormalization_fock(c_float k) const
    {
        assert(parent != nullptr);
        int index = parent->momentumRanges.momentum_to_floor_index(k);
        if(index >= parent->momentumRanges.size() - 1) {
            index = parent->momentumRanges.size() - 2;
        }
        return Utility::Numerics::linearly_interpolate(k, parent->momentumRanges[index], parent->momentumRanges[index + 1], 
                                renormalization_cache[index][0], renormalization_cache[index + 1][0]);
    }
    c_float PhononInteraction::renormalization_infinity(c_float k) const
    {
        assert(parent != nullptr);
        int index = parent->momentumRanges.momentum_to_floor_index(k);
        if(index >= parent->momentumRanges.size() - 1) {
            index = parent->momentumRanges.size() - 2;
        }
        return -Utility::Numerics::linearly_interpolate(k, parent->momentumRanges[index], parent->momentumRanges[index + 1], 
                                renormalization_cache[index][1], renormalization_cache[index + 1][1]);
    }
    c_float PhononInteraction::fock_correction(c_float k) const
    {
        assert(parent != nullptr);
        int index = parent->momentumRanges.momentum_to_floor_index(k);
        if(index >= parent->momentumRanges.size() - 1) {
            index = parent->momentumRanges.size() - 2;
        }
        return Utility::Numerics::linearly_interpolate(k, parent->momentumRanges[index], parent->momentumRanges[index + 1], 
                                renormalization_cache[index][2], renormalization_cache[index + 1][2]);
    }
}