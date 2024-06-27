#pragma once
#include "GlobalDefinitions.hpp"
#include "ModelAttributes.hpp"
#include "MomentumRanges.hpp"
#include <cmath>
#include <limits>
#include <map>
#include <boost/math/special_functions/pow.hpp>
#include <SymbolicOperators/WickTerm.hpp>
#include <Utility/InputFileReader.hpp>
#include <Utility/Numerics/Interpolation.hpp>
#include <boost/math/quadrature/gauss.hpp>
#include <Utility/ConstexprPower.hpp>

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
		static constexpr int n_interpolate = 4;
		static constexpr int shifted_index(int index) {
			return (index > n_interpolate / 2 - 1 ? index + 1 - n_interpolate / 2 : 0);
		}
	public:
		ModelAttributes<c_complex> Delta;
		std::vector<c_float> continuum_boundaries() const;

		c_float temperature{};
		c_float phonon_coupling{};
		c_float omega_debye{};
		const c_float fermi_energy{};

		static constexpr c_float CUT_OFF = std::numeric_limits<c_float>::epsilon();

		inline c_complex interpolate_delta(c_float k) const {
			const int index = momentumRanges.momentum_to_floor_index(k);
			if (index >= DISCRETIZATION - 1) // Assuming Delta(k) = 0 for k -> infinity
				return (index >= DISCRETIZATION ? c_complex{} : Delta[DISCRETIZATION - 1]);
			if (index < 0) // Assuming Delta(k) = 0 for k->0
				return c_complex{};
			return Utility::Numerics::interpolate_from_vector<n_interpolate>(k, momentumRanges, Delta, shifted_index(index));
			//return Utility::Numerics::interpolate_lagrange<2>(k, 
			//	std::array<c_complex, 2>{ momentumRanges.index_to_momentum(index), momentumRanges.index_to_momentum(index + 1) },
			//	std::array<c_complex, 2>{ Delta[index], Delta[index + 1] });
		};
		inline c_float interpolate_delta_n(c_float k) const {
			const int index = momentumRanges.momentum_to_floor_index(k);
			if (index >= DISCRETIZATION - 1) // Assuming Delta(k) = 0 for k -> infinity
				return (index >= DISCRETIZATION ? c_float{} : std::real(Delta[2 * DISCRETIZATION - 1]));
			if (index < 0) // Assuming Delta(k) = const for k->0
				return c_float{};
			return Utility::Numerics::interpolate_from_vector<n_interpolate>(k, momentumRanges, Delta, shifted_index(index), DISCRETIZATION);
			//return Utility::Numerics::interpolate_lagrange<2>(k, 
			//	std::array<c_complex, 2>{ momentumRanges.index_to_momentum(index), momentumRanges.index_to_momentum(index + 1) },
			//	std::array<c_complex, 2>{ std::real(Delta[index + DISCRETIZATION]), std::real(Delta[index + DISCRETIZATION + 1]) });
		};

		constexpr static c_float bare_dispersion(c_float k) {
			return 0.5 * k * k;
		};

		c_complex sc_expectation_value(c_float k) const;
		c_float occupation(c_float k) const;
		inline c_float delta_n(c_float k) const {
			if(is_zero(fermi_wavevector - k)) {
				return 0.5 - occupation(k);
			} 
			else if(fermi_wavevector - k > 0) {
				return 1. - occupation(k);
			}
			return -occupation(k);
		}

		c_float fock_energy(c_float k) const;

		inline c_float bare_dispersion_to_fermi_level(c_float k) const {
			return bare_dispersion(k) - fermi_energy;
		};

		inline c_float energy(c_float k) const {
			return sqrt(boost::math::pow<2>(dispersion_to_fermi_level(k)) + std::norm(interpolate_delta(k)));
		};

		inline c_float dispersion_to_fermi_level(c_float k) const
		{
			return bare_dispersion_to_fermi_level(k) + fock_energy(k) + interpolate_delta_n(k);
		};

		void iterationStep(const ParameterVector& initial_values, ParameterVector& result);
		c_float computeCoefficient(SymbolicOperators::Coefficient const& coeff, c_float first, c_float second) const;
		inline c_float computeCoefficient(SymbolicOperators::Coefficient const& coeff, c_float first) const
		{
			return computeCoefficient(coeff, first, fermi_wavevector);
		};
		std::map<SymbolicOperators::OperatorType, std::vector<c_complex>> get_expectation_values() const;

		c_float g_lower_bound(c_float k) const;
		c_float g_upper_bound(c_float k) const;

		inline std::string info() const {
			return "SCModel // [T g omega_D epsilon_F] = [" + std::to_string(temperature)
				+ " " + std::to_string(phonon_coupling) + " " + std::to_string(omega_debye)
				+ " " + std::to_string(fermi_energy) + "]";
		}

		c_float internal_energy() const;

		SCModel(ModelInitializer const& parameters);
		virtual ~SCModel() = default;

		c_float fermi_wavevector{};
		const c_float V_OVER_N{};

		MomentumRanges momentumRanges;

		static constexpr int n_gauss = 60;
#ifdef _screening
		template<class ExpectationValues>
		inline auto integral_screening(ExpectationValues const& expecs, c_float k) const
		{
			if(is_zero(k)){
				auto integrand = [&](c_float q){
					return expecs(q) * q * q / (_screening * _screening + q * q);
				};
				const c_float prefactor = 2. * PhysicalConstants::em_factor;

				return prefactor * momentumRanges.integrate(integrand);
			}

			auto integrand = [&](c_float q){
				const c_float k_diff{ q - k };
				const c_float k_sum{ q + k };
				return expecs(q) * q * std::log((_screening * _screening + k_sum * k_sum) / (_screening * _screening + k_diff * k_diff));
			};
			const c_float prefactor = 0.5 * PhysicalConstants::em_factor / k;

			return prefactor * momentumRanges.integrate(integrand);
			//const c_float k_f_offset{  };
			//return prefactor * 
			//	(
			//		boost::math::quadrature::gauss<double, n_gauss>::integrate(integrand, momentumRanges.K_MIN, fermi_wavevector - k_f_offset)
			//		//+ boost::math::quadrature::gauss<double, n_gauss>::integrate(integrand, fermi_wavevector - k_f_offset, fermi_wavevector + k_f_offset)
			//		+ boost::math::quadrature::gauss<double, n_gauss>::integrate(integrand, fermi_wavevector + k_f_offset, momentumRanges.K_MAX)
			//	);
		}
#endif

		static constexpr c_float OFF_SET_FACTOR = 0.03;
		template<class ExpectationValues>
		inline auto I_1(ExpectationValues const& expecs, c_float k) const {
			auto integrand = [&expecs, &k](c_float x) {
				const c_float tanh_x = std::tanh(0.5 * x);
				return expecs(k * tanh_x) * x * tanh_x / Utility::constexprPower<2>(std::cosh(0.5 * x));
				};

			const c_float alpha{ momentumRanges.K_MIN / k };//
			if(is_zero(1. - alpha)) return decltype(expecs(k)){};
			const c_float lower_bound = std::log( (1. + alpha) / (1. - alpha) );
			constexpr c_float upper_bound = 30; // a priori error estimation of machine epsilon

			const c_float prefactor = 0.5 * PhysicalConstants::em_factor * k;
			constexpr c_float CUT_OFF_CONSTANT = 1.16034524813597e-11;//1.73136903588109e-7;

			if(k > fermi_wavevector){
				const c_float middle_bound = 2. * std::atanh(fermi_wavevector / k);
				const c_float k_f_offset = OFF_SET_FACTOR * DISCRETIZATION  * momentumRanges.STEP;

				if(k > fermi_wavevector + k_f_offset){
					const c_float middle_middle_bound = 2. * std::atanh((fermi_wavevector + k_f_offset) / k);

					return prefactor * (CUT_OFF_CONSTANT * expecs(k)
						+ boost::math::quadrature::gauss<double, n_gauss>::integrate(integrand, lower_bound, middle_bound)
						+ boost::math::quadrature::gauss<double, n_gauss>::integrate(integrand, middle_bound, middle_middle_bound)
						+ boost::math::quadrature::gauss<double, n_gauss>::integrate(integrand, middle_middle_bound, upper_bound));
				}

				return prefactor * (CUT_OFF_CONSTANT * expecs(k)
					+ boost::math::quadrature::gauss<double, n_gauss>::integrate(integrand, lower_bound, middle_bound)
					+ boost::math::quadrature::gauss<double, n_gauss>::integrate(integrand, middle_bound, upper_bound));
			}

			return prefactor * (CUT_OFF_CONSTANT * expecs(k)
				+ boost::math::quadrature::gauss<double, n_gauss>::integrate(integrand, lower_bound, upper_bound));
		}

		template<class ExpectationValues>
		inline auto I_2(ExpectationValues const& expecs, c_float k) const {
			auto integrand = [&expecs, &k](c_float x) {
				const c_float coth_x = 1. / std::tanh(0.5 * x);
				return expecs(k * coth_x) * x * coth_x / Utility::constexprPower<2>(std::sinh(0.5 * x));
				};
			const c_float beta{ momentumRanges.K_MAX / k };
			if(is_zero(beta - 1.)) return decltype(expecs(k)){};
			const c_float lower_bound = std::log( (beta + 1.) / (beta - 1.) );
			constexpr c_float upper_bound = 30; // a priori error estimation of machine epsilon

			const c_float prefactor = 0.5 * PhysicalConstants::em_factor * k;
			constexpr c_float CUT_OFF_CONSTANT = 1.16034524813640e-11;//1.73136904981569e-7;

			if(k < fermi_wavevector){
				const c_float middle_bound = 2. * std::atanh(k / fermi_wavevector); // acoth(x) = atanh(1/x)
				const c_float k_f_offset = OFF_SET_FACTOR * DISCRETIZATION * momentumRanges.STEP;

				if(k < fermi_wavevector - k_f_offset){
					const c_float middle_middle_bound = 2. * std::atanh(k / (fermi_wavevector - k_f_offset));

					return prefactor * (CUT_OFF_CONSTANT * expecs(k)
						+ boost::math::quadrature::gauss<double, n_gauss>::integrate(integrand, lower_bound, middle_bound)
						+ boost::math::quadrature::gauss<double, n_gauss>::integrate(integrand, middle_bound, middle_middle_bound)
						+ boost::math::quadrature::gauss<double, n_gauss>::integrate(integrand, middle_middle_bound, upper_bound));
				}

				return prefactor * (CUT_OFF_CONSTANT * expecs(k)
					+ boost::math::quadrature::gauss<double, n_gauss>::integrate(integrand, lower_bound, middle_bound)
					+ boost::math::quadrature::gauss<double, n_gauss>::integrate(integrand, middle_bound, upper_bound));
			}

			return prefactor * (CUT_OFF_CONSTANT * expecs(k)
				+ boost::math::quadrature::gauss<double, n_gauss>::integrate(integrand, lower_bound, upper_bound));
		}

	private:
		c_float compute_fermiwavevector(c_float epsilon_F) const;

		inline c_float phonon_alpha(const c_float k) const {
			if(fermi_wavevector - k == 0.0) return k * k;
			if(is_zero(k)) return 2 * PhysicalConstants::em_factor * fermi_wavevector * fermi_wavevector;

			const c_float log_expr = std::log(std::abs((fermi_wavevector + k) / (fermi_wavevector - k)));
			const c_float k2 = k * k;
			const c_float kF2 = fermi_wavevector * fermi_wavevector;

			return k2 - PhysicalConstants::em_factor * (kF2 - k2) * log_expr / k;
		}

		inline auto phonon_beta(const c_float k, const c_float ALPHA) const {
			return phonon_alpha(k) + ALPHA;
		}
	};
}