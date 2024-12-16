#pragma once
#include "GlobalDefinitions.hpp"
#include "MomentumRanges.hpp"
#include <array>
#include <vector>
#include "SplineContainer.hpp"

namespace Continuum {
    class SCModel;

    class PhononInteraction {
    private:
        void compute_renormalization_table();
        void compute_singularities();
    
        std::array<c_float, 2> get_singularities(c_float k) const;
    public:
        void set_parent(SCModel const * const _parent);

        // Delta epsilon + omega_D (as described in the CUT)
		c_float alpha_CUT(c_float k, c_float k_prime) const;
		// Delta epsilon - omega_D (as described in the CUT)
		c_float beta_CUT(c_float k, c_float k_prime) const;

        c_float sc_channel_lower_bound(c_float k) const;
		c_float sc_channel_upper_bound(c_float k) const;

        c_float renormalization_fock(c_float k) const;
        c_float renormalization_infinity(c_float k) const; 
        c_float fock_correction(c_float k) const;

        /*template<class ExpectationValues>
         void update_fock_correction(ExpectationValues const & delta_n) { 
            for (MomentumIterator it(& momentumRanges); it < MomentumIterator::max_idx(); ++it) {
                renormalization_cache[it.idx][2] = fock_channel_integral(delta_n, it.k);
            } 
        } */
        template<class ExpectationValues>
		decltype(std::declval<ExpectationValues>()(c_float{})) sc_channel_integral(ExpectationValues const& expecs, c_float k) const;

        template<class ExpectationValues>
		decltype(std::declval<ExpectationValues>()(c_float{})) fock_channel_integral(ExpectationValues const& expecs, c_float k) const;
    private:
        SCModel const * parent;
        c_float const * ptr_rho_F;
        c_float const * ptr_fermi_wavevector;
        c_float const * ptr_omega_debye;
        c_float const * ptr_phonon_coupling;
        MomentumRanges const * ptr_momentumRanges;
        // [0] lower singularity; [1] upper singularity
        std::vector<std::array<c_float, 2>> singularities_cache;
    public:
        // [0] up to k_F; [1] up to K_MAX
        std::vector<std::array<c_float, 3>> renormalization_cache;
    };

    template<class ExpectationValues>
	decltype(std::declval<ExpectationValues>()(c_float{})) PhononInteraction::sc_channel_integral(ExpectationValues const& expecs, c_float k) const
    {
        assert(parent != nullptr);
		const c_float prefactor = (*ptr_phonon_coupling) / (2. * PI * PI * (*ptr_rho_F));
		auto integrand = [&expecs](c_float q) {
			return q * q * expecs(q);
			};
		return prefactor * ptr_momentumRanges->integrate(integrand, sc_channel_lower_bound(k), sc_channel_upper_bound(k));
	}

    template<class ExpectationValues>
	decltype(std::declval<ExpectationValues>()(c_float{})) PhononInteraction::fock_channel_integral(ExpectationValues const& expecs, c_float k) const
    {
#ifdef PHONON_SC_CHANNEL_ONLY
        return decltype(std::declval<ExpectationValues>()(c_float{})){};
#else
        assert(parent != nullptr);
        auto integrand_fock = [this, k, &expecs](c_float q) -> c_float {
	    	return expecs(k) * q * q * 
	    		( 1. / std::copysign(std::abs(beta_CUT(q, k)) + CUT_REGULARIZATION, beta_CUT(q, k)) 
	    		- 1. / std::copysign(std::abs(alpha_CUT(q,k)) + CUT_REGULARIZATION, alpha_CUT(q, k)) );
	    };
    
        const auto singularities = get_singularities(k);
        decltype(std::declval<ExpectationValues>()(c_float{})) result;
        if (singularities[0] < (*ptr_fermi_wavevector)) {
            result = ptr_momentumRanges->integrate(integrand_fock, ptr_momentumRanges->K_MIN, singularities[0]);

            if (singularities[1] < (*ptr_fermi_wavevector)) {
                result += ptr_momentumRanges->integrate(integrand_fock, singularities[0], singularities[1]) +
                        ptr_momentumRanges->integrate(integrand_fock, singularities[1], (*ptr_fermi_wavevector));
            }
            else {
                result += ptr_momentumRanges->integrate(integrand_fock, singularities[0], (*ptr_fermi_wavevector));
            }
        }
        else {
            result = ptr_momentumRanges->integrate(integrand_fock, ptr_momentumRanges->K_MIN, (*ptr_fermi_wavevector));
        }
        return result;
#endif
    }
}