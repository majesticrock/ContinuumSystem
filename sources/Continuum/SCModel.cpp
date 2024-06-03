#include "SCModel.hpp"
#include "../../../Utility/sources/Numerics/Integration/TrapezoidalRule.hpp"
#include "../../../Utility/sources/ConstexprPower.hpp"
#include "../../../Utility/sources/Numerics/Minimization/Bisection.hpp"
#include <algorithm>
#include <numeric>
#include <complex>

#include <boost/math/quadrature/gauss.hpp>

using Utility::constexprPower;

#define _use_coulomb

namespace Continuum {
#ifdef approximate_theta
	constexpr c_float delta_range_factor = 1;
	constexpr c_float sc_is_finite_range = 1;
#else
	constexpr c_float delta_range_factor = 15;
	constexpr c_float sc_is_finite_range = 15;
#endif

	SCModel::SCModel(ModelInitializer const& parameters)
		: Delta(DISCRETIZATION, parameters.phonon_coupling* parameters.omega_debye), temperature{ parameters.temperature },
		phonon_coupling{ parameters.phonon_coupling }, omega_debye{ parameters.omega_debye }, fermi_energy{ parameters.fermi_energy },
		fermi_wavevector{ sqrt(2 * parameters.fermi_energy) },
		V_OVER_N{ fermi_wavevector > 0 ? 3. * PI * PI / (constexprPower<3>(fermi_wavevector)) : 1 },
		K_MAX{ sqrt(2 * (fermi_energy + delta_range_factor * omega_debye)) - fermi_wavevector },
		K_MIN{ fermi_energy > delta_range_factor * omega_debye ? sqrt(2 * (fermi_energy - delta_range_factor * omega_debye)) - fermi_wavevector : -fermi_wavevector },
		STEP{ (K_MAX - K_MIN) / (DISCRETIZATION - 1) },
		MAX_K_WITH_SC{ sqrt(2 * (fermi_energy + sc_is_finite_range * omega_debye)) },
		MIN_K_WITH_SC{ sqrt(2 * (fermi_energy - sc_is_finite_range * omega_debye)) }
	{
		omega_debye += PRECISION;
		assert(index_to_momentum(0) >= 0);
		Delta = decltype(Delta)::FromAllocator([&](size_t i) -> c_complex {
			const c_float k = index_to_momentum(i);
			const c_float magnitude = (k < sqrt(2. * (fermi_energy - omega_debye)) || k > sqrt(2. * (fermi_energy + omega_debye)))
				? 0.01 : 0.1;
			if (i < DISCRETIZATION) {
				return std::polar(magnitude, i * 2.0 * PI / (DISCRETIZATION)+0.5 * PI);
			}
			return c_complex{};
			}, DISCRETIZATION);

		//std::cout << std::abs(bare_dispersion_to_fermi_level(index_to_momentum(0))) << std::endl;
		//Delta = decltype(Delta)::Random(2 * DISCRETIZATION);
		//Delta = decltype(Delta)::Gaussian(DISCRETIZATION, DISCRETIZATION / 2, (int)(DISCRETIZATION / delta_range_factor), 0.01 * phonon_coupling);
	}

	std::vector<c_float> SCModel::continuum_boundaries() const
	{
		return {
			2 * this->energy(Utility::Numerics::Minimization::bisection([this](c_float k) { return this->energy(k); },
				index_to_momentum(0), index_to_momentum(DISCRETIZATION - 1), PRECISION, 100)),
			2 * this->energy(Utility::Numerics::Minimization::bisection([this](c_float k) { return -this->energy(k); },
				index_to_momentum(0), index_to_momentum(DISCRETIZATION - 1), PRECISION, 100))
		};
	}

	c_float SCModel::dispersion_to_fermi_level(c_float k) const {
		return bare_dispersion_to_fermi_level(k);// + interpolate_delta_n(k); // TODO
	}

	c_complex SCModel::sc_expectation_value(c_float k) const {
		const auto DELTA = interpolate_delta(k);
		const c_float E = energy(k);
		if (is_zero(DELTA)) return 0;
		if (is_zero(temperature)) {
			return -DELTA / (2 * E);
		}
		return -std::tanh(E / (2 * temperature)) * DELTA / (2 * E);
	}

	c_float SCModel::occupation(c_float k) const {
		const auto DELTA = std::abs(interpolate_delta(k));
		const auto eps_mu = dispersion_to_fermi_level(k);
		if (is_zero(DELTA)) {
			if (is_zero(temperature)) {
				if (is_zero(eps_mu)) return 0.5;
				return (eps_mu < 0 ? 1 : 0);
			}
			return 1. / (1 + std::exp(eps_mu / temperature));
		}
		const c_float E = sqrt(eps_mu * eps_mu + DELTA * DELTA);
		if (is_zero(temperature)) {
			return 0.5 * (1 - (eps_mu / E));
		}
		return 0.5 * (1 - (eps_mu / E) * std::tanh(E / (2 * temperature)));
	}

	void SCModel::iterationStep(const ParameterVector& initial_values, ParameterVector& result) {
		result.setZero();
		this->Delta.fill_with(initial_values);
		auto phonon_integrand = [this](c_float x) -> c_complex {
			return x * x * sc_expectation_value(x);
			};

#ifdef _use_coulomb
		auto sc_em_inner_integrand = [this](c_float k_tilde) -> c_complex {
			return k_tilde * sc_expectation_value(k_tilde);
			};
		/* auto num_em_inner_integrand = [this](c_float k_tilde) -> c_float {
			return k_tilde * occupation(k_tilde);
			}; */
#endif

			//#pragma omp parallel for
		for (int u_idx = 0; u_idx < DISCRETIZATION; ++u_idx) {
			const c_float k = index_to_momentum(u_idx);
#ifdef approximate_theta
			// approximate theta(omega - 0.5*|l^2 - k^2|) as theta(omega - eps_k)theta(omega - eps_l)
			if (std::abs(bare_dispersion_to_fermi_level(k)) > omega_debye) {
				continue;
			}
#endif

#ifdef _use_coulomb
			auto sc_em_integrand = [&](c_float q) -> c_complex {
				//if(is_zero(q)) {
				//	return 2 * k * sc_expectation_value(k);
				//}
				const c_float lower = std::max(std::abs(q - k), MIN_K_WITH_SC);
				const c_float upper = std::min(q + k, MAX_K_WITH_SC);

				if (lower >= upper) return c_complex{};
				return boost::math::quadrature::gauss<double, 30>::integrate(
					sc_em_inner_integrand, lower, upper) * (q / (q * q + fermi_wavevector * fermi_wavevector));
				};

			/*auto num_em_integrand = [&](c_float q) -> c_float {
				if(is_zero(q)) {
					return k * occupation(k);
				}
				const c_float lower = std::abs(q - k);
				const c_float upper = q + k;

				const c_float numerics_lower = std::max(lower, MIN_K_WITH_SC);
				const c_float numerics_upper = std::min(upper, MAX_K_WITH_SC);

				c_float result{};
				// if Delta = 0 => <n_k> = 1 if epsilon < 0 and 0 otherwise
				if(lower < numerics_lower) {
					result = 0.5 * (numerics_lower * numerics_lower - lower * lower);
				}
				if(numerics_lower < numerics_upper)
					result += boost::math::quadrature::gauss<double, 30>::integrate(
						num_em_inner_integrand, numerics_lower, numerics_upper) / q;

				return result;
			};*/

			result(u_idx) = (PhysicalConstants::em_factor / k) * boost::math::quadrature::gauss<double, 30>::integrate(
				sc_em_integrand, c_float{}, MAX_K_WITH_SC + k);
			//result(u_idx + DISCRETIZATION) = (-PhysicalConstants::em_factor / k) * boost::math::quadrature::gauss<double, 30>::integrate(
			//		num_em_integrand, c_float{}, MAX_K_WITH_SC + k);
#endif
			result(u_idx) -= (V_OVER_N * phonon_coupling / (2. * PI * PI))
				* boost::math::quadrature::gauss<double, 30>::integrate(
					phonon_integrand, g_lower_bound(k), g_upper_bound(k));
		}
		this->Delta.fill_with(result);
		result -= initial_values;
	}

	c_float SCModel::computeCoefficient(SymbolicOperators::Coefficient const& coeff, c_float first, c_float second) const
	{
		if (coeff.name == "\\epsilon_0")
		{
			return bare_dispersion_to_fermi_level(first);
		}
		else if (coeff.name == "g")
		{
#ifdef approximate_theta
			if (omega_debye > std::abs(bare_dispersion_to_fermi_level(first))
				&& omega_debye > std::abs(bare_dispersion_to_fermi_level(second)))
#else
			if (coeff.momenta[0] == coeff.momenta[1])
			{
				return this->phonon_coupling * this->V_OVER_N;
			}
			if (omega_debye > std::abs(bare_dispersion(first) - bare_dispersion(second)))
#endif
			{
				return this->phonon_coupling * this->V_OVER_N;
			}
			else
			{
				return c_float{};
			}
		}
		else
		{
			throw std::invalid_argument("Coefficient not recognized! " + coeff.name);
		}
	}

	std::map<SymbolicOperators::OperatorType, std::vector<c_complex>> SCModel::get_expectation_values() const
	{
		std::map<SymbolicOperators::OperatorType, std::vector<c_complex>> ret;
		ret.emplace(SymbolicOperators::Number_Type, std::vector<c_complex>(DISCRETIZATION));
		ret.emplace(SymbolicOperators::SC_Type, std::vector<c_complex>(DISCRETIZATION));

		for (int k = 0; k < DISCRETIZATION; ++k) {
			const c_float momentum = index_to_momentum(k);
			ret.at(SymbolicOperators::Number_Type)[k] = this->occupation(momentum);
			ret.at(SymbolicOperators::SC_Type)[k] = this->sc_expectation_value(momentum);
		}

		return ret;
	}

	c_float SCModel::internal_energy() const {
		auto procedure = [this](c_float k) -> c_float {
			return -k * k * energy(k);
			};
		return (V_OVER_N / (2 * PI * PI)) * boost::math::quadrature::gauss<c_float, 30>::integrate(
			procedure, MIN_K_WITH_SC, MAX_K_WITH_SC);
	}
}