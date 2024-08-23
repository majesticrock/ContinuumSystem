#include "MomentumRanges.hpp"
#include <iostream>

#ifdef approximate_theta
constexpr Continuum::c_float inner_offset = 1. - 1e-5;
#else
#ifndef mielke_coulomb
constexpr Continuum::c_float inner_offset = 10;
#else
constexpr Continuum::c_float inner_offset = 2;
#endif
#endif

namespace Continuum {
	MomentumRanges::MomentumRanges(c_float* k_F, const c_float omega_debye)
		: 
#ifndef mielke_coulomb
		K_MAX{ 2 * (*k_F) }, K_MIN{ 0 },
#else
		K_MAX{ (*k_F) + (inner_offset * 2) * omega_debye }, K_MIN{ (*k_F) - (inner_offset * 2) * omega_debye },
#endif
#ifndef approximate_theta

		INNER_K_MAX{ (*k_F) + inner_offset * omega_debye / (*k_F) },
		INNER_K_MIN{ (*k_F) - inner_offset * omega_debye / (*k_F) },
#else
		INNER_K_MAX{ sqrt((*k_F**k_F) + 2 * inner_offset * omega_debye) },
		INNER_K_MIN{ sqrt((*k_F**k_F) - 2 * inner_offset * omega_debye) },
#endif
		LOWER_STEP{ (INNER_K_MIN - K_MIN) / _OUTER_DISC },
		INNER_STEP{ (INNER_K_MAX - INNER_K_MIN) / _INNER_DISC },
		UPPER_STEP{ (K_MAX - INNER_K_MAX) / _OUTER_DISC },
		K_F{ k_F }
	{ 
		assert(index_to_momentum(0) >= 0);
#ifndef _NDEBUG
		for(int i = 0; i < DISCRETIZATION; ++i) {
			const double k = index_to_momentum(i);
			if(i != momentum_to_index(k)){
				std::cerr << "MomentumRanges.cpp:" << __LINE__ 
					<< " || " << i << " -> " << k << " -> " << momentum_to_index(k) << std::endl;
				throw -1;
			} 
		}
#endif
	}

	c_float MomentumRanges::index_to_momentum(int k_idx) const
	{
		if(k_idx < _OUTER_DISC) {
			return K_MIN + LOWER_STEP * k_idx;
		}
		k_idx -= _OUTER_DISC;
		if(k_idx < _INNER_DISC) {
			return INNER_K_MIN + INNER_STEP * k_idx;
		}
		k_idx -= _INNER_DISC;
		return INNER_K_MAX + UPPER_STEP * k_idx;
	}

	int MomentumRanges::momentum_to_index(c_float k) const
	{
		if(k < INNER_K_MIN) {
			return static_cast<int>(std::lround((k - K_MIN) / LOWER_STEP));
		}
		if(k < INNER_K_MAX) {
			return static_cast<int>(std::lround((k - INNER_K_MIN) / INNER_STEP) + _OUTER_DISC);
		}
		return static_cast<int>(std::lround((k - INNER_K_MAX) / UPPER_STEP) + _INNER_DISC + _OUTER_DISC);
	}

	int MomentumRanges::momentum_to_floor_index(c_float k) const
	{
		if(k < INNER_K_MIN) {
			return static_cast<int>((k - K_MIN) / LOWER_STEP);
		}
		if(k < INNER_K_MAX) {
			return static_cast<int>((k - INNER_K_MIN) / INNER_STEP + _OUTER_DISC);
		}
		return static_cast<int>((k - INNER_K_MAX) / UPPER_STEP + _INNER_DISC + _OUTER_DISC);
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