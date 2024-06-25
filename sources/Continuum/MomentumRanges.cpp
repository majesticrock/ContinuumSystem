#include "MomentumRanges.hpp"
#include <iostream>

constexpr Continuum::c_float sc_is_finite_range = 1;

//#define _uniform_disc

namespace Continuum {
    MomentumRanges::MomentumRanges(const c_float k_F, const c_float omega_debye)
        : K_MAX{ 2 * k_F }, // + sc_is_finite_range * sqrt(2 * omega_debye)
		K_MIN{ 0 }, //k_F - sc_is_finite_range * sqrt(2 * omega_debye)
		INNER_K_MAX{ k_F + 5 * omega_debye },
		INNER_K_MIN{ k_F - 5 * omega_debye },
        STEP{ 
#ifdef _uniform_disc
            (K_MAX - K_MIN) / DISCRETIZATION
#else
            (INNER_K_MIN - K_MIN) / _OUTER_DISC
#endif
        },
		INNER_STEP{ (INNER_K_MAX - INNER_K_MIN) / _INNER_DISC },
        K_F{ k_F }
    { 
        assert(index_to_momentum(0) >= 0);

        for(int i = 0; i < DISCRETIZATION; ++i) {
            const double k = index_to_momentum(i);
            if(i != momentum_to_index(k)){
                std::cerr << i << " -> " << k << " -> " << momentum_to_index(k) << std::endl;
            } 
        }
    }

    c_float MomentumRanges::index_to_momentum(int k_idx) const
    {
#ifdef _uniform_disc
        return K_MIN + STEP * k_idx;
#else
        if(k_idx < _OUTER_DISC) {
            return K_MIN + STEP * k_idx;
        }
        k_idx -= _OUTER_DISC;
        if(k_idx < _INNER_DISC) {
            return INNER_K_MIN + INNER_STEP * k_idx;
        }
        k_idx -= _INNER_DISC;
        return INNER_K_MAX + STEP * k_idx;
#endif
    }

    int MomentumRanges::momentum_to_index(c_float k) const
    {
#ifdef _uniform_disc
        return static_cast<int>(std::lround((k - K_MIN) / STEP));
#else
        if(k < INNER_K_MIN) {
            return static_cast<int>(std::lround((k - K_MIN) / STEP));
        }
        if(k < INNER_K_MAX) {
            return static_cast<int>(std::lround((k - INNER_K_MIN) / INNER_STEP) + _OUTER_DISC);
        }
        return static_cast<int>(std::lround((k - INNER_K_MAX) / STEP) + _INNER_DISC + _OUTER_DISC);
#endif
    }

    std::vector<c_float> MomentumRanges::get_k_points() const
    {
        std::vector<c_float> ks;
		ks.resize(DISCRETIZATION);
		for (int i = 0; i < DISCRETIZATION; ++i) {
			ks[i] = index_to_momentum(i);
		}
		return ks;
    }
}