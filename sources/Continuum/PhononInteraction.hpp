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

        c_float renormalization_flow(c_float k) const;
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
        std::vector<std::array<c_float, 2>> renormalization_cache;
    };


    /*
    *  Template implementations
    */
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

        const c_float M_SQUARED = (*ptr_phonon_coupling) * (*ptr_omega_debye) / (*ptr_rho_F);
        const c_float factor = M_SQUARED / (2. * PI * PI);
        const auto singularities = get_singularities(k);

        auto integrand_fock_alpha = [this, &k, &expecs](c_float q) -> c_float {
            return expecs(q) * q * q / std::copysign(std::abs(alpha_CUT(q, k)) + CUT_REGULARIZATION, alpha_CUT(q, k));
		};
        auto integrand_fock_beta = [this, &k, &expecs](c_float q) -> c_float {
            return expecs(q) * q * q / std::copysign(std::abs(beta_CUT(q, k)) + CUT_REGULARIZATION, beta_CUT(q, k));
		};
        
        decltype(std::declval<ExpectationValues>()(c_float{})) result = factor *
            ( ptr_momentumRanges->cpv_integrate(integrand_fock_alpha, singularities[0])
            + ptr_momentumRanges->cpv_integrate(integrand_fock_beta,  singularities[1]) );
        return result;
#endif
    }
}