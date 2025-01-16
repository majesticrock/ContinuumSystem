#include "PhononInteraction.hpp"
#include "SCModel.hpp"
#include <mrock/utility/Numerics/Roots/Bisection.hpp>
#include <mrock/utility/Numerics/Interpolation.hpp>
#include <cmath>

namespace Continuum {
    void PhononInteraction::compute_renormalization_table()
    {
        assert(parent != nullptr);
        renormalization_cache.resize(parent->momentumRanges.size(), std::array<c_float, 2>{ c_float{}, c_float{} });

#ifdef NO_FOCK_PHONON
        return;
#else
        // Typically around 1e-3
        const c_float M_SQUARED = (*ptr_phonon_coupling) * (*ptr_omega_debye) / (*ptr_rho_F);
        const c_float factor = M_SQUARED / (2. * PI * PI);
#ifdef CUT_DISPERSION_NO_COULOMB
        // TODO: Stimmt das Vorzeichen?!
        for (MomentumIterator it(& parent->momentumRanges); it < MomentumIterator::max_idx(); ++it) {
            // Fock correction is initially 0
            const c_float sqrt_minus = sqrt(std::abs(it.k * it.k - *ptr_omega_debye));
            const c_float sqrt_plus = sqrt(it.k * it.k + *ptr_omega_debye);

            const c_float term_alpha = (it.k * it.k < *ptr_omega_debye 
                ? sqrt_minus * (std::atan(*ptr_fermi_wavevector / sqrt_minus) - M_PI_2)
                : sqrt_minus * std::log( (std::abs(*ptr_fermi_wavevector - sqrt_minus) + CUT_REGULARIZATION) / (std::abs(*ptr_fermi_wavevector + sqrt_minus) + CUT_REGULARIZATION) ) );
            const c_float term_beta = sqrt_plus * std::log( (std::abs(*ptr_fermi_wavevector - sqrt_plus) + CUT_REGULARIZATION) / (std::abs(*ptr_fermi_wavevector + sqrt_plus) + CUT_REGULARIZATION) );
            renormalization_cache[it.idx][0] = factor * ( term_alpha - term_beta );

        }
#else
        auto compute_at_k = [this](c_float k) {
            // First singularity is in alpha, second one in beta
            const std::array<c_float, 2> singularities = get_singularities(k);
            auto integrand_fock_alpha = [this, &k](c_float q) -> c_float {
                // for q->infinity the integrand becomes q^2 / (0.5 * q^2) = 2
                // to avoid giant values we substract this 2
                // This is essentially just a constant shift in energy, which we simply absorb into the chemical potential
                return c_float{-2} + q * q / std::copysign(std::abs(alpha_CUT(q, k)) + CUT_REGULARIZATION, alpha_CUT(q, k));
		    };
            auto integrand_fock_beta = [this, &k](c_float q) -> c_float {
                return c_float{-2} + q * q / std::copysign(std::abs(beta_CUT(q, k)) + CUT_REGULARIZATION, beta_CUT(q, k));
		    };

            return mrock::utility::Numerics::Integration::GeneralizedPrincipalValue<c_float, 120>::generalized_principal_value(
                    integrand_fock_alpha, parent->fermi_wavevector, parent->momentumRanges.K_MAX /* infinity */, singularities[0])
                + boost::math::quadrature::gauss<c_float, 10000>::integrate(integrand_fock_alpha, parent->momentumRanges.K_MAX, 1e4 /* infinity */)
                + mrock::utility::Numerics::Integration::GeneralizedPrincipalValue<c_float, 120>::generalized_principal_value(
                    integrand_fock_beta, parent->momentumRanges.K_MIN, parent->fermi_wavevector, singularities[1]);
        };

        const c_float value_at_k_f = compute_at_k(parent->fermi_wavevector);
        for (MomentumIterator it(& parent->momentumRanges); it < MomentumIterator::max_idx(); ++it) {
            // Fock correction is initially 0
            renormalization_cache[it.idx][0] = factor * ( compute_at_k(it.k) - value_at_k_f );
        }
#endif
#endif
    }

    void PhononInteraction::compute_singularities()
    {
        assert(parent != nullptr);
        singularities_cache.resize(parent->momentumRanges.size(), std::array<c_float, 2>{ c_float{}, c_float{} });
        for (MomentumIterator it(& parent->momentumRanges); it < MomentumIterator::max_idx(); ++it) 
        {
            try {
			    singularities_cache[it.idx][0] = mrock::utility::Numerics::Roots::bisection(
                    [this, &it](c_float q) { return alpha_CUT(q, it.k);  }, parent->momentumRanges.K_MIN, parent->momentumRanges.K_MAX, PRECISION, 200);
            } catch (mrock::utility::Numerics::Roots::NoRootException const & e) {
                singularities_cache[it.idx][0] = 2 * parent->momentumRanges.K_MAX;
            }
            try {
                singularities_cache[it.idx][1] = mrock::utility::Numerics::Roots::bisection(
                    [this, &it](c_float q) { return beta_CUT(q, it.k); }, parent->momentumRanges.K_MIN, parent->momentumRanges.K_MAX, PRECISION, 200);
            } catch (mrock::utility::Numerics::Roots::NoRootException const & e) {
                singularities_cache[it.idx][1] = 2 * parent->momentumRanges.K_MAX;
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
            mrock::utility::Numerics::linearly_interpolate(k, parent->momentumRanges[index], parent->momentumRanges[index + 1], 
                                singularities_cache[index][0], singularities_cache[index + 1][0]),
            mrock::utility::Numerics::linearly_interpolate(k, parent->momentumRanges[index], parent->momentumRanges[index + 1], 
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
		return mrock::utility::Numerics::Roots::bisection(func, parent->momentumRanges.K_MIN, fermi_wavevector, PRECISION, 200);
#else
		const auto lb = func(parent->momentumRanges.K_MIN);
		const auto ub = func(k);
		if (lb * ub >= c_float{}) return parent->momentumRanges.K_MIN;
		return mrock::utility::Numerics::Roots::bisection(func, parent->momentumRanges.K_MIN, k, PRECISION, 200);
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
		return mrock::utility::Numerics::Roots::bisection(func, fermi_wavevector, parent->momentumRanges.K_MAX, PRECISION, 200);
#else
		const auto lb = func(k);
		const auto ub = func(parent->momentumRanges.K_MAX);
		if (lb * ub >= c_float{}) return parent->momentumRanges.K_MAX;
		return mrock::utility::Numerics::Roots::bisection(func, k, parent->momentumRanges.K_MAX, PRECISION, 200);
#endif
	}

    c_float PhononInteraction::renormalization_flow(c_float k) const
    {
        assert(parent != nullptr);
        int index = parent->momentumRanges.momentum_to_floor_index(k);
        if(index >= parent->momentumRanges.size() - 1) {
            index = parent->momentumRanges.size() - 2;
        }
        return mrock::utility::Numerics::linearly_interpolate(k, parent->momentumRanges[index], parent->momentumRanges[index + 1], 
                                renormalization_cache[index][0], renormalization_cache[index + 1][0]);
    }

    c_float PhononInteraction::fock_correction(c_float k) const
    {
        assert(parent != nullptr);
        int index = parent->momentumRanges.momentum_to_floor_index(k);
        if(index >= parent->momentumRanges.size() - 1) {
            index = parent->momentumRanges.size() - 2;
        }
        return mrock::utility::Numerics::linearly_interpolate(k, parent->momentumRanges[index], parent->momentumRanges[index + 1], 
                                renormalization_cache[index][1], renormalization_cache[index + 1][1]);
    }

    c_float PhononInteraction::fock_channel(c_float k, c_float k_prime) const
    {
        assert(parent != nullptr);
        // TODO: Think about the factor of 2
        const c_float factor = (*ptr_phonon_coupling) * (*ptr_omega_debye) / (*ptr_rho_F);
        return factor * ( std::copysign(std::abs(beta_CUT(k, k_prime))  + CUT_REGULARIZATION, beta_CUT(k, k_prime))
                        - std::copysign(std::abs(alpha_CUT(k, k_prime)) + CUT_REGULARIZATION, alpha_CUT(k, k_prime)) );
    }
}