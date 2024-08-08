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

		inline c_float bare_dispersion(c_float k) const;
		inline c_float bare_dispersion_to_fermi_level(c_float k) const;
		inline c_float dispersion_to_fermi_level(c_float k) const;
		inline c_float dispersion_to_fermi_level(int k) const;
		inline c_float energy(c_float k) const;
		inline c_float energy(int k) const;
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
		std::vector<c_float> continuum_boundaries() const;
		std::vector<c_complex> phonon_gap() const;
		std::vector<c_complex> coulomb_gap() const;
		std::vector<c_float> single_particle_dispersion() const;
		const std::map<SymbolicOperators::OperatorType, std::vector<c_complex>>& get_expectation_values() const;

		void set_new_parameters(ModelInitializer const& parameters);
		SCModel(ModelInitializer const& parameters);
		virtual ~SCModel() = default;

		template<class ExpectationValues>
		decltype(std::declval<ExpectationValues>()(c_float{})) integral_phonon(ExpectationValues const& expecs, c_float k) const
		{
			auto integrand = [&](c_float q) {
				return q * q * expecs(q);
				};
			return -(phonon_coupling / (2. * PI * PI)) * fermi_wavevector * fermi_wavevector * fermi_wavevector
				* boost::math::quadrature::gauss<double, 60>::integrate(integrand, g_lower_bound(k), g_upper_bound(k));
		}

		template<class ExpectationValues>
		decltype(std::declval<ExpectationValues>()(c_float{})) integral_screening(ExpectationValues const& expecs, c_float k) const
		{
#ifndef mielke_coulomb
			if (is_zero(k)) {
				auto integrand = [&](c_float q) {
					return expecs(q) * q * q / (_screening * _screening + q * q);
					};
				const c_float prefactor = 2. * PhysicalConstants::em_factor * coulomb_scaling * fermi_wavevector;
				return prefactor * momentumRanges.integrate(integrand);
			}
			auto integrand = [&](c_float q) {
				const c_float k_diff{ q - k };
				const c_float k_sum{ q + k };
				return expecs(q) * q * std::log((_screening * _screening + k_sum * k_sum) / (_screening * _screening + k_diff * k_diff));
				};
			const c_float prefactor = 0.5 * PhysicalConstants::em_factor * coulomb_scaling * fermi_wavevector / k;
			return prefactor * momentumRanges.integrate(integrand);
#else
			auto integrand = [&](c_float q) {
				return expecs(q) * q * q;
				};
			return coulomb_scaling * PhysicalConstants::em_factor * momentumRanges.integrate(integrand);
#endif
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

		void set_splines();
	public:
		ModelAttributes<c_complex> Delta;
		c_float temperature{};
		c_float phonon_coupling{};
		c_float omega_debye{};
		c_float coulomb_scaling{};
		c_float r_s{};

		c_float fermi_wavevector{};
		c_float fermi_energy{};

		SplineContainer occupation;
		SplineContainer sc_expectation_value;

		MomentumRanges momentumRanges;
	private:
		mutable std::map<SymbolicOperators::OperatorType, std::vector<c_complex>> _expecs;
	};

	// /////////// //
	//   Inlines   //
	// /////////// //
	c_float SCModel::bare_dispersion(c_float k) const {
		return (0.5 / PhysicalConstants::effective_mass) * fermi_wavevector * fermi_wavevector * k * k;
	}
	c_float SCModel::delta_n(c_float k) const {
		if (is_zero(1. - k)) {
			return 0.5 - std::real(occupation(k));
		}
		else if (1. - k > 0) {
			return 1. - std::real(occupation(k));
		}
		return -std::real(occupation(k));
	}
	c_float SCModel::delta_n_index(int k) const {
		if (is_zero(1. - momentumRanges.index_to_momentum(k))) {
			return 0.5 - occupation_index(k);
		}
		else if (1. - momentumRanges.index_to_momentum(k) > 0) {
			return 1. - occupation_index(k);
		}
		return -occupation_index(k);
	}
	c_float SCModel::bare_dispersion_to_fermi_level(c_float k) const {
		return bare_dispersion(k) - fermi_energy;
	}
	c_float SCModel::dispersion_to_fermi_level(c_float k) const {
		return bare_dispersion_to_fermi_level(k) + fock_energy(k) + interpolate_delta_n(k);
	}
	c_float SCModel::dispersion_to_fermi_level(int k) const {
		return bare_dispersion_to_fermi_level(momentumRanges.index_to_momentum(k)) + fock_energy(momentumRanges.index_to_momentum(k)) + std::real(Delta[k + DISCRETIZATION]);
	}
	c_float SCModel::energy(c_float k) const {
		return sqrt(boost::math::pow<2>(dispersion_to_fermi_level(k)) + std::norm(interpolate_delta(k)));
	}
	c_float SCModel::energy(int k) const {
		return sqrt(boost::math::pow<2>(dispersion_to_fermi_level(k)) + std::norm(Delta[k]));
	}
}