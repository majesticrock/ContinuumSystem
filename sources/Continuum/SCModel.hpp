#pragma once
#include "GlobalDefinitions.hpp"
#include "ModelAttributes.hpp"
#include <cmath>
#include <limits>
#include <boost/math/special_functions/pow.hpp>

namespace Continuum {
	struct ModelInitializer {
		c_float temperature;
		c_float U;
		c_float omega_debye;
	};

	class SCModel {
	public:
		ModelAttributes<c_complex> Delta;

		c_float temperature{};
		c_float U{};
		c_float omega_debye{};

		static constexpr c_float CUT_OFF = std::numeric_limits<c_float>::epsilon();

		inline c_complex interpolate_delta(c_float k) const {
			const int index = momentum_to_index(k);
			if (index + 1 == DISCRETIZATION)
				return Delta[index]; // Assuming Delta(k) = const for k -> infinity
			const c_float k_lower = k - index_to_momentum(index);
			const c_float k_upper = index_to_momentum(index + 1) - k;
			// k_lower + k_upper = k_n+1 - k_n
			return (Delta[index] * k_upper + Delta[index + 1] * k_lower) / (k_upper + k_lower);
		};

		inline c_float energy(c_float k) const {
			return sqrt(boost::math::pow<4>(k) + std::norm(interpolate_delta(k)));
		};
		inline c_complex sc_expectation_value(c_float k) const {
			const c_float E = energy(k);
			if (is_zero(E)) return 0;
			if (is_zero(temperature)) {
				return interpolate_delta(k) / (2 * E);
			}
			return std::tanh(E / (2 * temperature)) * interpolate_delta(k) / (2 * E);
		};
		void iterationStep(const ParameterVector& initial_values, ParameterVector& result);

		SCModel(ModelInitializer const& parameters);
	};
}