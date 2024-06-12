#pragma once
#include "GlobalDefinitions.hpp"
#include "ModelAttributes.hpp"
#include <cmath>
#include <limits>
#include <map>
#include <boost/math/special_functions/pow.hpp>
#include "../../../FermionCommute/sources/WickTerm.hpp"
#include <Utility/InputFileReader.hpp>
#include <Utility/Numerics/Interpolation.hpp>
#include <boost/math/quadrature/gauss.hpp>

namespace Continuum {
	struct ModelInitializer {
		c_float temperature;
		c_float phonon_coupling;
		c_float omega_debye;
		c_float fermi_energy;

		ModelInitializer(Utility::InputFileReader& input)
			: temperature{ PhysicalConstants::k_B * input.getDouble("T") }, phonon_coupling{ input.getDouble("phonon_coupling") },
			omega_debye{ input.getDouble("omega_debye") }, fermi_energy{ input.getDouble("fermi_energy") }
		{ };
	};

	class SCModel {
	public:
		ModelAttributes<c_complex> Delta;

		inline c_float index_to_momentum(int k_idx) const {
			assert(fermi_wavevector + K_MIN + STEP * k_idx >= 0);
			return fermi_wavevector + K_MIN + STEP * k_idx;
		}
		inline int momentum_to_index(c_float k) const {
			assert(k >= 0);
			return static_cast<int>(std::lround((k - fermi_wavevector - K_MIN) / STEP));
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
		c_float phonon_coupling{};
		c_float omega_debye{};
		const c_float fermi_energy{};

		static constexpr c_float CUT_OFF = std::numeric_limits<c_float>::epsilon();

		inline c_complex interpolate_delta(c_float k) const {
			const int index = momentum_to_index(k);
			if (index >= DISCRETIZATION - 1) // Assuming Delta(k) = 0 for k -> infinity
				return (index >= DISCRETIZATION ? c_complex{} : Delta[DISCRETIZATION - 1]);
			if (index < 0) // Assuming Delta(k) = 0 for k->0
				return c_complex{};
			return Utility::Numerics::linearly_interpolate(k, index_to_momentum(index), index_to_momentum(index + 1),
				Delta[index], Delta[index + 1]);
		};
		inline c_float interpolate_delta_n(c_float k) const {
			const int index = momentum_to_index(k);
			if (index >= DISCRETIZATION - 1) // Assuming Delta(k) = 0 for k -> infinity
				return (index >= DISCRETIZATION ? c_float{} : std::real(Delta[2 * DISCRETIZATION - 1]));
			if (index < 0) // Assuming Delta(k) = const for k->0
				return c_float{};
			return Utility::Numerics::linearly_interpolate(k, index_to_momentum(index), index_to_momentum(index + 1),
				std::real(Delta[index + DISCRETIZATION]), std::real(Delta[index + DISCRETIZATION + 1]));
		};
		constexpr static c_float bare_dispersion(c_float k) {
			return 0.5 * k * k;
		};

	public:
		c_complex sc_expectation_value(c_float k) const;
		c_float occupation(c_float k) const;

		inline c_float fock_energy(k) const {
			if(is_zero(k - fermi_wavevector)) {
				return -PhysicalConstants::em_factor * fermi_wavevector;
			}

			return -PhysicalConstants::em_factor * fermi_wavevector * (
				1.0 + ((fermi_wavevector * fermi_wavevector - k * k) / (2.0 * k * fermi_wavevector)) 
					* std::log(std::abs((k + fermi_wavevector) / (k - fermi_wavevector)))
			);
		}

		inline c_float bare_dispersion_to_fermi_level(c_float k) const {
			return bare_dispersion(k) - fermi_energy;
		};
		c_float dispersion_to_fermi_level(c_float k) const;
		inline c_float energy(c_float k) const {
			return sqrt(boost::math::pow<2>(dispersion_to_fermi_level(k)) + std::norm(interpolate_delta(k)));
		};

		void iterationStep(const ParameterVector& initial_values, ParameterVector& result);
		c_float computeCoefficient(SymbolicOperators::Coefficient const& coeff, c_float first, c_float second) const;
		inline c_float computeCoefficient(SymbolicOperators::Coefficient const& coeff, c_float first) const
		{
			return computeCoefficient(coeff, first, fermi_wavevector);
		};
		std::map<SymbolicOperators::OperatorType, std::vector<c_complex>> get_expectation_values() const;

		inline c_float g_lower_bound(c_float k) const {
#ifdef approximate_theta
			return sqrt(std::max((static_cast<c_float>(2) * (fermi_energy - omega_debye)), c_float{}));
#else
			return sqrt(std::max(c_float{}, k * k - 2 * omega_debye));
#endif
		}
		inline c_float g_upper_bound(c_float k) const {
#ifdef approximate_theta
			return sqrt(static_cast<c_float>(2) * (fermi_energy + omega_debye));
#else
			return sqrt(k * k + 2 * omega_debye);
#endif
		}
		inline std::string info() const {
			return "SCModel // [T g omega_D epsilon_F] = [" + std::to_string(temperature)
				+ " " + std::to_string(phonon_coupling) + " " + std::to_string(omega_debye)
				+ " " + std::to_string(fermi_energy) + "]";
		}

		c_float internal_energy() const;

		SCModel(ModelInitializer const& parameters);
		virtual ~SCModel() = default;

		const c_float fermi_wavevector{};
		const c_float V_OVER_N{};

		const c_float K_MAX{};
		const c_float K_MIN{};
		const c_float STEP{};

		const c_float MAX_K_WITH_SC;
		const c_float MIN_K_WITH_SC;

	private:
		template<class ExpectationValues>
		inline auto I_1(ExpecationValues const& expecs, c_float k) const {
			auto integrand = [&expecs, &k](c_float y) {
				if(is_zero(y)) return decltype(expecs(k)){};
				return (1. - 2. * y) * std::log(1. / y - 1.) * expecs(k * (1. - 2. * y));
				};
			constexpr c_float lower_bound = 0.0;
			const c_float upper_bound = 0.5 * (1. - MIN_K_WITH_SC / k);
			if(is_zero(lower_bound - upper_bound)) return decltype(expecs(k)){};

			return 2. * PhysicalConstants::em_prefactor * k * boost::math::quadrature::gauss<double, 30>::integrate(integrand, lower_bound, upper_bound);
		}
		template<class ExpectationValues>
		inline auto I_2(ExpecationValues const& expecs, c_float k) const {
			auto integrand = [&expecs, &k](c_float y) {
				if(is_zero(y)) return decltype(expecs(k)){};
				return (2. * y + 1.) * std::log(1. + 1. / y) * expecs(k * (2. * y + 1.));
				};
			constexpr c_float lower_bound = 0.0;
			const c_float upper_bound = 0.5 * (MAX_K_WITH_SC / k - 1.0);
			if(is_zero(lower_bound - upper_bound)) return decltype(expecs(k)){};

			return 2. * PhysicalConstants::em_prefactor * k * boost::math::quadrature::gauss<double, 30>::integrate(integrand, lower_bound, upper_bound);
		}
	};
}