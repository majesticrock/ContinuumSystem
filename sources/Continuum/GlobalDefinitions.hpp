#pragma once
#define _USE_MATH_DEFINES
#include <Eigen/Dense>
#include <complex>
#include <cmath>
#include <type_traits>
#include "../../../Utility/sources/UnderlyingFloatingPoint.hpp"

namespace Continuum {
	typedef double c_float;
	typedef std::complex<c_float> c_complex;

	using SpinorMatrix = Eigen::Matrix<c_complex, Eigen::Dynamic, Eigen::Dynamic>;
	using ParameterVector = Eigen::Vector<c_complex, Eigen::Dynamic>;

	constexpr c_complex I = { 0, 1 };

	template<class RealType>
	constexpr RealType PRECISION = 1e2 * std::numeric_limits<RealType>::epsilon();
	template<> constexpr float PRECISION<float> = 1e1f * std::numeric_limits<float>::epsilon();

	template<class NumberType>
	constexpr bool is_zero(const NumberType& number) {
		return abs(number) < PRECISION<Utility::UnderlyingFloatingPoint_t<NumberType>>;
	}

	constexpr c_float U_MAX = 1;
	constexpr int DISCRETIZATION = 500;
	constexpr c_float STEP = U_MAX / DISCRETIZATION;

	constexpr c_float index_to_momentum(int u) {
		return STEP * u / (1 - STEP * u);
	}
	inline int momentum_to_index(c_float k) {
		return static_cast<int>(k / (STEP * (k + 1)));
	}
}