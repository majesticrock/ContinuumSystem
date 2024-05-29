#pragma once
#define _USE_MATH_DEFINES
#define BOOST_MATH_GAUSS_NO_COMPUTE_ON_DEMAND

#include <Eigen/Dense>
#include <complex>
#include <cmath>
#include <type_traits>
#include <cstdint>
#include <cstring>
#include <limits>
#include "../../../Utility/sources/UnderlyingFloatingPoint.hpp"

//#define approximate_theta

namespace Continuum {
	using c_float = double;
	using c_complex = std::complex<c_float>; //

	constexpr c_float PI = static_cast<c_float>(M_PIl);

	using SpinorMatrix = Eigen::Matrix<c_complex, Eigen::Dynamic, Eigen::Dynamic>;
	using ParameterVector = Eigen::Vector<c_complex, Eigen::Dynamic>;

	//constexpr c_complex I = { 0, 1 };

	/* We use 
	*  hbar = 1
	*  m_e = 1
	*  e = 1
	*/
	namespace PhysicalConstants {
		constexpr c_float k_B = 8.617333262e-5; // eV / K
		constexpr c_float vacuum_permitivity = 0.05526349406 * 3.62262628; // sqrt(eV)
		constexpr c_float em_factor = 1. / (2 * PI * PI * vacuum_permitivity); // sqrt(eV)
	}

	constexpr double SQRT_PRECISION = 9.5367431640625e-07; 
	constexpr double PRECISION = 9.0949470177292824e-13; // 0 | 01111101011 | 0000000000000000000000000000000000000000000000000000

	/* This function abuses the structure of our desired precision:
	*  The mantissa is empty, i.e., we can solely rely on the exponent.
	*  If the exponent of <number> is >= 0b01111101011, |number| >= precision, i.e., not 0
	*/
	inline bool is_zero(double number) {
		static_assert( std::numeric_limits<double>::is_iec559, "IEEE 754 floating point not verified!" );
		// 9.0949470177292824e-13   <->    0 | 01111101011 | 0000000000000000000000000000000000000000000000000000
		uint64_t tmp; // silence compiler warnings (at 0 overhead)
		std::memcpy(&tmp, &number, sizeof(tmp));
		return static_cast<uint16_t>((tmp >> 52) & 0x7FF) < 0b01111010111;
	}

	inline bool is_zero(std::complex<double> number) {
		return is_zero(std::abs(number));
	}

	extern int DISCRETIZATION;
	extern c_float INV_N;

	inline void set_discretization(int N) {
		DISCRETIZATION = N;
		INV_N = 1.0 / N;
	}
}