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
		c_float fermi_energy;

		ModelInitializer(Utility::InputFileReader& input)
			: temperature{ PhysicalConstants::k_B * input.getDouble("T") }, U{ input.getDouble("U") },
			omega_debye{ input.getDouble("omega_debye") }, fermi_energy{ input.getDouble("fermi_energy") }
		{ };
	};

	class SCModel {
	public:
		ModelAttributes<c_complex> Delta;

		inline c_float index_to_momentum(int u) const {
			assert(fermi_wavevector + U_MIN + STEP * u >= 0);
			return fermi_wavevector + U_MIN + STEP * u;
		}
		inline int momentum_to_index(c_float k) const {
			assert(k >= 0);
			return static_cast<int>(std::lround((k - fermi_wavevector - U_MIN) / STEP));
		}
		inline std::vector<c_float> get_k_points() const {
			std::vector<c_float> ks;
			ks.resize(DISCRETIZATION);
			for (int i = 0; i < DISCRETIZATION; ++i) {
				ks[i] = index_to_momentum(i);
			}
			return ks;
		}
		std::vector<c_float> continuum_boundaries() const;
	protected:
		c_float temperature{};
		c_float U{};
		c_float omega_debye{};
		c_float fermi_energy{};
		c_float fermi_wavevector{};

		c_float V_OVER_N{};

		static constexpr c_float CUT_OFF = std::numeric_limits<c_float>::epsilon();

		inline c_complex interpolate_delta(c_float k) const {
			const int index = momentum_to_index(k);
			if (index >= DISCRETIZATION - 1)
				return Delta.back(); // Assuming Delta(k) = const for k -> infinity
			if (index < 0)
				return Delta.front();
			return Utility::Numerics::linearly_interpolate(k, index_to_momentum(index), index_to_momentum(index + 1), Delta[index], Delta[index + 1]);
		};

		constexpr static c_float bare_dispersion(c_float k) {
			return 0.5 * k * k;
		};
		inline c_float bare_dispersion_to_fermi_level(c_float k) const {
			return bare_dispersion(k) - fermi_energy;
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
#ifdef approximate_theta
			return sqrt(std::max((2. * (fermi_energy - omega_debye)), c_float{}));
#else
			return sqrt(std::max(c_float{}, k * k - 2 * omega_debye));
#endif
		}
		inline c_float u_upper_bound(c_float k) const {
#ifdef approximate_theta
			return sqrt(2. * (fermi_energy + omega_debye));
#else
			return sqrt(k * k + 2 * omega_debye);
#endif
		}
		inline std::string info() const {
			return "SCModel // [T U omega_D mu] = [" + std::to_string(temperature)
				+ " " + std::to_string(U) + " " + std::to_string(omega_debye)
				+ " " + std::to_string(fermi_energy) + "]";
		}

		SCModel(ModelInitializer const& parameters);

		c_float U_MAX{};
		c_float U_MIN{};
		c_float STEP{};
	};
}