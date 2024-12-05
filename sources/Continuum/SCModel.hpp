#pragma once
#include "GlobalDefinitions.hpp"
#include "ModelAttributes.hpp"
#include "MomentumRanges.hpp"
#include "SplineContainer.hpp"
#include "ModelInitializer.hpp"
#include <cmath>
#include <limits>
#include <utility>
#include <map>
#include <boost/math/special_functions/pow.hpp>
#include <SymbolicOperators/WickTerm.hpp>
#include <Utility/InputFileReader.hpp>
#include <Utility/Numerics/Interpolation.hpp>
#include <boost/math/quadrature/gauss.hpp>
#include <Utility/ConstexprPower.hpp>

namespace Continuum {
	class SCModel {
	public:
		c_complex interpolate_delta(c_float k) const;
		c_float interpolate_delta_n(c_float k) const;

		inline c_float bare_dispersion_to_fermi_level(c_float k) const;
		inline c_float dispersion_to_fermi_level(c_float k) const;
		inline c_float dispersion_to_fermi_level_index(int k) const;
		inline c_float energy(c_float k) const;
		inline c_float energy_index(int k) const;
		c_float fock_energy(c_float k) const;

		inline c_float delta_n(c_float k) const;

		c_float g_lower_bound(c_float k) const;
		c_float g_upper_bound(c_float k) const;

		c_complex k_infinity_integral() const;
		c_complex k_zero_integral() const;

		void iterationStep(const ParameterVector& initial_values, ParameterVector& result);
		inline c_float computeCoefficient(SymbolicOperators::Coefficient const& coeff, c_float first) const {
			return computeCoefficient(coeff, first, fermi_wavevector);
		}
		c_float computeCoefficient(SymbolicOperators::Coefficient const& coeff, c_float first, c_float second) const;

		std::string info() const;
		std::string to_folder() const;
		c_float internal_energy() const;
		std::vector<c_complex> phonon_gap() const;
		std::vector<c_complex> coulomb_gap() const;
		c_complex delta_max() const;
		std::vector<c_float> single_particle_dispersion() const;
		const std::map<SymbolicOperators::OperatorType, std::vector<c_complex>>& get_expectation_values() const;

		void set_new_parameters(ModelInitializer const& parameters);
		SCModel(ModelInitializer const& parameters);
		virtual ~SCModel() = default;

		template<class ExpectationValues>
		decltype(std::declval<ExpectationValues>()(c_float{})) integral_phonon(ExpectationValues const& expecs, c_float k) const
		{
			const c_float prefactor = phonon_coupling / (2. * PI * PI * rho_F);
			auto integrand = [&expecs](c_float q) {
				return q * q * expecs(q);
				};

			return prefactor * momentumRanges.integrate(integrand, g_lower_bound(k), g_upper_bound(k));
		}

		template<class ExpectationValues>
		decltype(std::declval<ExpectationValues>()(c_float{})) integral_screening(ExpectationValues const& expecs, c_float k) const
		{
			if (is_zero(coulomb_scaling)) return decltype(std::declval<ExpectationValues>()(c_float{})) {};
			if (is_zero(k)) {
				auto integrand = [&](c_float q) {
					return expecs(q) * q * q / (screening * screening + q * q);
					};
				const c_float prefactor = 2. * coulomb_scaling * PhysicalConstants::em_factor;
				return prefactor * momentumRanges.integrate(integrand);
			}
			auto integrand = [&](c_float q) {
				const c_float k_diff{ q - k };
				const c_float k_sum{ q + k };
				return expecs(q) * q * std::log((screening * screening + k_sum * k_sum) / (screening * screening + k_diff * k_diff));
				};
			const c_float prefactor = 0.5 * coulomb_scaling * PhysicalConstants::em_factor / k;
			return prefactor * momentumRanges.integrate(integrand);
		}

	private:
		static constexpr int n_interpolate = 4;
		static constexpr int shifted_index(int index) {
			return (index > n_interpolate / 2 - 1 ? index + 1 - n_interpolate / 2 : 0);
		}

		c_complex sc_expectation_value_index(int k) const;
		c_float occupation_index(int k) const;
		inline c_float delta_n_index(int k) const;

		// Effectively 2 * ( bare_dispersion + E_fock )
		c_float phonon_alpha(const c_float k) const;
		inline auto phonon_beta(const c_float k, const c_float ALPHA) const {
			return phonon_alpha(k) - ALPHA;
		}

		inline c_float log_expression(c_float k_sum, c_float k_diff) const {
			return std::log((screening * screening + k_sum * k_sum) / (screening * screening + k_diff * k_diff));
		}

		// Calls the publically available fock_energy() if we consider the full Coulomb interaction and returns 0 otherwise
		inline c_float __fock_energy(c_float k) const {
#ifdef _coulomb_only_in_bcs_channel
			return c_float{};
#else
			return this->fock_energy(k);
#endif
		}

		// Calls the publically available interpolate_delta_n() if we consider the full Coulomb interaction and returns 0 otherwise
		inline c_float __interpolate_delta_n(c_float k) const {
#ifdef _coulomb_only_in_bcs_channel
			return c_float{};
#else
			return this->interpolate_delta_n(k);
#endif
		}

		void set_splines();
	public:
		ModelAttributes<c_complex> Delta;
		c_float temperature{};
		c_float phonon_coupling{};
		c_float omega_debye{};
		c_float fermi_energy{};
		c_float screening_ratio{};
		c_float screening{};
		c_float coulomb_scaling{};

		SplineContainer occupation;
		SplineContainer sc_expectation_value;

		c_float fermi_wavevector{};
		c_float rho_F{};

		MomentumRanges momentumRanges;
	private:
		mutable std::map<SymbolicOperators::OperatorType, std::vector<c_complex>> _expecs;
	};

	// /////////// //
	//   Inlines   //
	// /////////// //
	c_float SCModel::delta_n(c_float k) const {
		if (is_zero(fermi_wavevector - k)) {
			return 0.5 - std::real(occupation(k));
		}
		else if (fermi_wavevector - k > 0) {
			return 1. - std::real(occupation(k));
		}
		return -std::real(occupation(k));
	}
	c_float SCModel::delta_n_index(int k) const {
		if (is_zero(fermi_wavevector - momentumRanges.index_to_momentum(k))) {
			return 0.5 - occupation_index(k);
		}
		else if (fermi_wavevector - momentumRanges.index_to_momentum(k) > 0) {
			return 1. - occupation_index(k);
		}
		return -occupation_index(k);
	}
	c_float SCModel::bare_dispersion_to_fermi_level(c_float k) const {
		return bare_dispersion(k) - fermi_energy;
	}
	c_float SCModel::dispersion_to_fermi_level(c_float k) const {
		return bare_dispersion_to_fermi_level(k) + __fock_energy(k) + __interpolate_delta_n(k);
	}
	c_float SCModel::dispersion_to_fermi_level_index(int k) const {
#ifdef _coulomb_only_in_bcs_channel
		return bare_dispersion_to_fermi_level(momentumRanges.index_to_momentum(k));
#else
		return bare_dispersion_to_fermi_level(momentumRanges.index_to_momentum(k)) + __fock_energy(momentumRanges.index_to_momentum(k)) + std::real(Delta[k + DISCRETIZATION]);
#endif
	}
	c_float SCModel::energy(c_float k) const {
		return sqrt(boost::math::pow<2>(dispersion_to_fermi_level(k)) + std::norm(interpolate_delta(k)));
	}
	c_float SCModel::energy_index(int k) const {
		return sqrt(boost::math::pow<2>(dispersion_to_fermi_level_index(k)) + std::norm(Delta[k]));
	}
}