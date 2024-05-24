#include "SCModel.hpp"
#include "../../../Utility/sources/Numerics/Integration/TrapezoidalRule.hpp"
#include "../../../Utility/sources/ConstexprPower.hpp"
#include "../../../Utility/sources/Numerics/Minimization/Bisection.hpp"
#include <algorithm>
#include <numeric>

#include <boost/math/quadrature/gauss_kronrod.hpp>

using Utility::constexprPower;

namespace Continuum {
	SCModel::SCModel(ModelInitializer const& parameters)
		: Delta(DISCRETIZATION, parameters.phonon_coupling* parameters.omega_debye), temperature{ parameters.temperature },
		phonon_coupling{ parameters.phonon_coupling }, omega_debye{ parameters.omega_debye }, fermi_energy{ parameters.fermi_energy },
		fermi_wavevector{ sqrt(2 * parameters.fermi_energy) },
		V_OVER_N{ fermi_wavevector > 0 ? 3. * PI * PI / (constexprPower<3>(fermi_wavevector)) : 1 },
		K_MAX{ sqrt(2 * (fermi_energy + omega_debye)) - fermi_wavevector },
		K_MIN{ fermi_energy > omega_debye ? sqrt(2 * (fermi_energy - omega_debye)) - fermi_wavevector : -fermi_wavevector },
		STEP{ (K_MAX - K_MIN) / (DISCRETIZATION - 1) }
	{
		omega_debye += SQRT_PRECISION<c_float>;
		assert(index_to_momentum(0) >= 0);
		//std::cout << std::abs(bare_dispersion_to_fermi_level(index_to_momentum(0))) << std::endl;
		//Delta = decltype(Delta)::Random(DISCRETIZATION);
	}

	std::vector<c_float> SCModel::continuum_boundaries() const
	{
		return {
			2 * this->energy(Utility::Numerics::Minimization::bisection([this](c_float k) { return this->energy(k); },
				index_to_momentum(0), index_to_momentum(DISCRETIZATION - 1), PRECISION<c_float>, 100)),
			2 * this->energy(Utility::Numerics::Minimization::bisection([this](c_float k) { return -this->energy(k); },
				index_to_momentum(0), index_to_momentum(DISCRETIZATION - 1), PRECISION<c_float>, 100))
		};
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
		const auto eps_mu = bare_dispersion_to_fermi_level(k);
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

	c_complex SCModel::sc_from_epsilon(c_float epsilon) const
	{
		const auto DELTA = std::abs(interpolate_delta(sqrt(2.0 * (epsilon + fermi_energy))));
		const c_float E = sqrt(epsilon * epsilon + DELTA * DELTA);//energy(k);
		if (is_zero(DELTA)) return 0;
		if (is_zero(temperature)) {
			return -DELTA / (2 * E);
		}
		return -std::tanh(E / (2 * temperature)) * DELTA / (2 * E);
	}

	c_float SCModel::n_from_epsilon(c_float epsilon) const
	{
		const auto DELTA = std::abs(interpolate_delta(sqrt(2.0 * (epsilon + fermi_energy))));
		if (is_zero(DELTA)) {
			if (is_zero(temperature)) {
				if (is_zero(epsilon)) return 0.5;
				return (epsilon < 0 ? 1 : 0);
			}
			return 1. / (1 + std::exp(epsilon / temperature));
		}
		const c_float E = sqrt(epsilon * epsilon + DELTA * DELTA);
		if (is_zero(temperature)) {
			return 0.5 * (1 - (epsilon / E));
		}
		return 0.5 * (1 - (epsilon / E) * std::tanh(E / (2 * temperature)));
	}

	void SCModel::iterationStep(const ParameterVector& initial_values, ParameterVector& result) {
		result.setZero();
		this->Delta.fill_with(initial_values);

#ifndef _use_epsilon_integration
		auto integrand = [this](c_float x) -> c_complex {
			return x * x * sc_expectation_value(x);
			};

		//result.fill(Utility::Numerics::Integration::trapezoidal_rule_kahan(integrand, g_lower_bound(0), g_upper_bound(0), DISCRETIZATION));
		c_float error;
		result.fill(boost::math::quadrature::gauss_kronrod<c_float, 61>::integrate(integrand, g_lower_bound(0), g_upper_bound(0), 8, 1e-14, &error));
#else
		auto integrand = [this](c_float eps) -> c_complex {
			return sqrt(2 * (eps + this->fermi_energy)) * sc_from_epsilon(eps);
			};
		c_float error;
		result.fill(boost::math::quadrature::gauss_kronrod<c_float, 61>::integrate(integrand, -omega_debye, omega_debye, 8, 1e-14, &error));
#endif

//		for (int u_idx = 0; u_idx < DISCRETIZATION; ++u_idx) {
//			const c_float k = index_to_momentum(u_idx);
//#ifdef approximate_theta
//			// approximate theta(omega - 0.5*|l^2 - k^2|) as theta(omega - eps_k)theta(omega - eps_l)
//			if (std::abs(bare_dispersion_to_fermi_level(k)) > omega_debye) {
//				continue;
//			}
//#endif
//			double error;
//			result(u_idx) = boost::math::quadrature::gauss_kronrod<double, 61>::integrate(integrand, g_lower_bound(k), g_upper_bound(k), 8, 1e-16, &error);
//			//result(u_idx) = Utility::Numerics::Integration::trapezoidal_rule(integrand, g_lower_bound(k), g_upper_bound(k), 1000);
//		}

		const c_float prefactor = -V_OVER_N / (2. * PI * PI);
		result *= prefactor * phonon_coupling;
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
}