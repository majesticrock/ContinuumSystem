#pragma once
#include "GlobalDefinitions.hpp"
#include "ModelAttributes.hpp"
#include <cmath>
#include <limits>
#include <map>
#include <boost/math/special_functions/pow.hpp>
#include "../../../FermionCommute/sources/WickTerm.hpp"
#include "../../../Utility/sources/InputFileReader.hpp"
#include "../../../Utility/sources/Numerics/Interpolation.hpp"

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
	protected:
		c_float temperature{};
		c_float U{};
		c_float omega_debye{};
		c_float chemical_potential{};

		static constexpr c_float CUT_OFF = std::numeric_limits<c_float>::epsilon();

		inline c_complex interpolate_delta(c_float k) const {
			const int index = momentum_to_index(k);
			if (index + 1 >= DISCRETIZATION)
				return Delta[index]; // Assuming Delta(k) = const for k -> infinity
			return Utility::Numerics::linearly_interpolate(k, index_to_momentum(index), index_to_momentum(index + 1), Delta[index], Delta[index + 1]);
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

	public:
		void iterationStep(const ParameterVector& initial_values, ParameterVector& result);
		c_float computeCoefficient(SymbolicOperators::Coefficient const& coeff, c_float first, c_float second = c_float{}) const;
		std::map<SymbolicOperators::OperatorType, std::vector<c_complex>> get_expectation_values() const;

		inline c_float u_lower_bound(c_float k) const {
			return sqrt(std::max(c_float{}, k * k - 2 * omega_debye));
		}
		inline c_float u_upper_bound(c_float k) const {
			return sqrt(k * k + 2 * omega_debye);
		}
		inline std::string info() const {
			return "SCModel // [T U omega_D mu] = [" + std::to_string(temperature)
				+ " " + std::to_string(U) + " " + std::to_string(omega_debye)
				+ " " + std::to_string(chemical_potential) + "]";
		}

		SCModel(ModelInitializer const& parameters);
	};
}