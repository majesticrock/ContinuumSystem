#pragma once
#include "GlobalDefinitions.hpp"
#include "ModelAttributes.hpp"
#include <cmath>
#include <limits>
#include <boost/math/special_functions/pow.hpp>
#include "../../../FermionCommute/sources/WickTerm.hpp"
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
	protected:
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
		c_complex sc_expectation_value(c_float k) const;
		c_float occupation(c_float k) const;
		void iterationStep(const ParameterVector& initial_values, ParameterVector& result);

	public:
		inline std::string info() const {
			return "SCModel // [T U omega_D mu] = [" + std::to_string(temperature)
				+ " " + std::to_string(U) + " " + std::to_string(omega_debye)
				+ " " + std::to_string(chemical_potential) + "]";
		}

		c_float computeCoefficient(SymbolicOperators::WickTerm const& term, c_float k, c_float l) const;
		c_complex computeTerm(SymbolicOperators::WickTerm const& term, c_float k, c_float l) const;

		SCModel(ModelInitializer const& parameters);
	};
}