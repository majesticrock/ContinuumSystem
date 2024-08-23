#pragma once
#define _USE_MATH_DEFINES
//#define BOOST_MATH_GAUSS_NO_COMPUTE_ON_DEMAND

#include <Eigen/Dense>
#include <complex>
#include <cmath>
#include <type_traits>
#include <cstdint>
#include <cstring>
#include <limits>
#include <Utility/UnderlyingFloatingPoint.hpp>
#include <iostream>

//#define approximate_theta
//#define mielke_coulomb
//#define _complex
#define _screening 1e-4

namespace Continuum {
	using c_float = double;
#ifdef _complex
	using c_complex = std::complex<c_float>;
#else
	using c_complex = c_float;
#endif

	constexpr c_float PI = static_cast<c_float>(M_PI);
	constexpr c_float PI_2 = static_cast<c_float>(M_PI_2);

	constexpr c_float SQRT_PRECISION = 2.384185791015625e-07;
	constexpr c_float PRECISION = 5.684341886080802e-14; // 0 | 01111010011 | 0000000000000000000000000000000000000000000000000000

	/* This function abuses the structure of our desired precision:
	*  The mantissa is empty, i.e., we can solely rely on the exponent.
	*  If the exponent of <number> is >= 0b01111010011, |number| >= precision, i.e., not 0
	*/
	inline bool is_zero(double number) {
		static_assert(std::numeric_limits<double>::is_iec559, "IEEE 754 floating point not verified!");
		// 5.684341886080802e-14   <->   0 | 01111010011 | 0000000000000000000000000000000000000000000000000000
		uint64_t tmp; // silence compiler warnings (at 0 overhead)
		std::memcpy(&tmp, &number, sizeof(tmp));
		return static_cast<uint16_t>((tmp >> 52) & 0x7FF) < 0b01111010011;
	}

	inline bool is_zero(std::complex<double> number) {
		return is_zero(std::abs(number));
	}

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
		constexpr c_float vacuum_permitivity = 0.05526349406 * 3.62262628; // 0.2001989859063799115 sqrt(eV)
		constexpr c_float screening_prefactor = 0.4107320221286488672; // eV^(1/4) -> sqrt(1 / (3 * pi * pi * epsilon_0))
		constexpr c_float em_factor = 1. / (4 * PI * PI * vacuum_permitivity); // 1 / (4 * pi * pi * epsilon_0) = 0.12652559550141668 sqrt(eV)
		constexpr c_float effective_mass = 1; // m* / m_e - lead: 2.1
		// 0.1265255955014166767 sqrt(eV)
	}

	constexpr c_float bare_dispersion(c_float k) {
		return (0.5 / PhysicalConstants::effective_mass) * k * k;
	};
	
	extern int DISCRETIZATION;
	extern c_float INV_N;

	static constexpr int REL_INNER_DISCRETIZATION = 2;
	extern int _INNER_DISC;
	extern int _OUTER_DISC;

	inline void set_discretization(int N) {
		DISCRETIZATION = N;
		INV_N = 1.0 / N;
		_INNER_DISC = N / REL_INNER_DISCRETIZATION;
		_OUTER_DISC = (N - _INNER_DISC) / 2;
	}
}