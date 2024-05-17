#pragma once
#define _USE_MATH_DEFINES
#include <Eigen/Dense>
#include <complex>
#include <cmath>
#include <type_traits>
#include "../../../Utility/sources/UnderlyingFloatingPoint.hpp"

#define approximate_theta

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
	constexpr RealType SQRT_PRECISION = 5e-6; //5e1 * 1.4901174450357931e-8; // 50 * sqrt(epsilon)
	template<> constexpr float SQRT_PRECISION<float> = 1.09182874e-3f; // sqrt(10 * epsilon)

	template<class RealType>
	constexpr RealType PRECISION = SQRT_PRECISION<RealType> *SQRT_PRECISION<RealType>;

	template<class NumberType>
	constexpr bool is_zero(const NumberType& number) {
		return abs(number) < PRECISION<Utility::UnderlyingFloatingPoint_t<NumberType>>;
	}

	extern int DISCRETIZATION;
	extern c_float INV_N;

	inline void set_discretization(int N) {
		DISCRETIZATION = N;
		INV_N = 1.0 / N;
	}
}