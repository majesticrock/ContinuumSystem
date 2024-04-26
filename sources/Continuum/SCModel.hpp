#pragma once
#include "GlobalDefinitions.hpp"
#include "ModelAttributes.hpp"
#include <cmath>
#include <limits>
#include <boost/math/special_functions/pow.hpp>

#include "../../../Utility/sources/InputFileReader.hpp"

namespace Continuum {
	struct ModelInitializer {
		c_float temperature;
		c_float U;
		c_float omega_debye;
		c_float chemical_potential;

		ModelInitializer(Utility::InputFileReader& input)
			: temperature{ PhysicalConstants::k_B * input.getDouble("T") }, U{ input.getDouble("U") },
			omega_debye{ input.getDouble("omega_debye") }, chemical_potential{ input.getDouble("chemical_potential") }
		{ };
	};

	class SCModel {
	public:
		ModelAttributes<c_complex> Delta;

		c_float temperature{};
		c_float U{};
		c_float omega_debye{};
		c_float chemical_potential{};

		static constexpr c_float CUT_OFF = std::numeric_limits<c_float>::epsilon();

		inline c_complex interpolate_delta(c_float k) const {
			const int index = momentum_to_index(k);
			if (index + 1 >= DISCRETIZATION)
				return Delta[index]; // Assuming Delta(k) = const for k -> infinity
			const c_float k_upper = index_to_momentum(index + 1) - k;
			const c_float k_lower = k - index_to_momentum(index);
			// k_lower + k_upper = k_n+1 - k_n
			return (Delta[index] * k_upper + Delta[index + 1] * k_lower) / (k_upper + k_lower);
		};
		constexpr static c_float bare_dispersion(c_float k) {
			return 0.5 * k * k;
		};
		inline c_float bare_dispersion_to_fermi_level(c_float k) const {
			return bare_dispersion(k) - chemical_potential;
		};
		inline c_float energy(c_float k) const {
			return sqrt(boost::math::pow<2>(bare_dispersion_to_fermi_level(k)) + std::norm(interpolate_delta(k)));
		};
		inline c_complex sc_expectation_value(c_float k) const {
			const auto DELTA = interpolate_delta(k);
			const c_float E = energy(k);
			if (is_zero(DELTA)) return 0;
			if (is_zero(temperature)) {
				return -DELTA / (2 * E);
			}
			return -std::tanh(E / (2 * temperature)) * DELTA / (2 * E);
		};
		inline c_float occupation(c_float k) const {
			const auto DELTA = interpolate_delta(k);
			if (is_zero(DELTA)) {
				if (is_zero(temperature)) {
					return (bare_dispersion_to_fermi_level(k) < 0 ? 1 : 0);
				}
				return 1. / (1 + std::exp(bare_dispersion_to_fermi_level(k) / temperature));
			}
			const c_float E = energy(k);
			if (is_zero(temperature)) {
				return bare_dispersion_to_fermi_level(k) / (2 * E);
			}
			return bare_dispersion_to_fermi_level(k) * std::tanh(E / (2 * temperature)) / (2 * E);
		};
		void iterationStep(const ParameterVector& initial_values, ParameterVector& result);
		inline std::string info() const {
			return "SCModel // [T U omega_D mu] = [" + std::to_string(temperature)
				+ " " + std::to_string(U) + " " + std::to_string(omega_debye)
				+ " " + std::to_string(chemical_potential) + "]";
		}

		SCModel(ModelInitializer const& parameters);
	};
}