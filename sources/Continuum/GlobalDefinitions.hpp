#pragma once
#define _USE_MATH_DEFINES
#include <Eigen/Dense>
#include <complex>
#include <cmath>
#include <type_traits>
#include "../../../Utility/sources/UnderlyingFloatingPoint.hpp"

#define _equal_disc

namespace Continuum {
	using c_float = double;
	using c_complex = std::complex<c_float>;

	using SpinorMatrix = Eigen::Matrix<c_complex, Eigen::Dynamic, Eigen::Dynamic>;
	using ParameterVector = Eigen::Vector<c_complex, Eigen::Dynamic>;

	constexpr c_complex I = { 0, 1 };

	namespace PhysicalConstants {
		constexpr c_float k_B = 8.617333262e-2; // meV / K
		constexpr c_float hbar = 4.135667696e-12; // meV s
		constexpr c_float electron_mass = 510998950; // meV / c^2
	}

	template<class RealType>
	constexpr RealType PRECISION = 1e2 * std::numeric_limits<RealType>::epsilon();
	template<> constexpr float PRECISION<float> = 1e1f * std::numeric_limits<float>::epsilon();

	template<class NumberType>
	constexpr bool is_zero(const NumberType& number) {
		return abs(number) < PRECISION<Utility::UnderlyingFloatingPoint_t<NumberType>>;
	}

	constexpr c_float U_MAX = 
#ifdef _equal_disc
	250;
#else
	1;
#endif
	constexpr int DISCRETIZATION = 1000;
	constexpr c_float STEP = U_MAX / DISCRETIZATION;

	constexpr c_float index_to_momentum(int u) {
#ifdef _equal_disc
		return STEP * u;
#else
		return STEP * u / (1 - STEP * u);
#endif
	}
	inline int momentum_to_index(c_float k) {
#ifdef _equal_disc
		return static_cast<int>(k / STEP);
#else
		return static_cast<int>(k / (STEP * (k + 1)));
#endif
	}
	inline std::vector<c_float> get_k_points() {
		std::vector<c_float> ks;
		ks.resize(DISCRETIZATION);
		for(int i = 0; i < DISCRETIZATION; ++i){
			ks[i] = index_to_momentum(i);
		}
		return ks;
	}
}