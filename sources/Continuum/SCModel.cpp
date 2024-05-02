#include "SCModel.hpp"
#include "../../../Utility/sources/Numerics/TrapezoidalRule.hpp"
#include <algorithm>

//#define approximate_theta

namespace Continuum {
	SCModel::SCModel(ModelInitializer const& parameters)
		: Delta(DISCRETIZATION, parameters.U), temperature{ parameters.temperature }, U{ parameters.U },
		omega_debye{ parameters.omega_debye }, chemical_potential{ parameters.chemical_potential }
	{ }

	c_complex sc_expectation_value(c_float k) const {
		const auto DELTA = interpolate_delta(k);
		const c_float E = energy(k);
		if (is_zero(DELTA)) return 0;
		if (is_zero(temperature)) {
			return -DELTA / (2 * E);
		}
		return -std::tanh(E / (2 * temperature)) * DELTA / (2 * E);
	}

	c_float occupation(c_float k) const {
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
	}

	void SCModel::iterationStep(const ParameterVector& initial_values, ParameterVector& result) {
		result.setZero();
		this->Delta.fill_with(initial_values);

		constexpr int N_L = 1000;

		for (int u_idx = 0; u_idx < DISCRETIZATION; ++u_idx) {
			const c_float k = index_to_momentum(u_idx);
#ifdef approximate_theta
			// approximate theta(omega - 0.5*|l^2 - k^2|) as theta(omega - 0.5*l^2)theta(omega - 0.5*k^2)
			if (bare_dispersion(k) - chemical_potential > omega_debye)
				continue;
			const c_float lower_bound = 0;
			const c_float upper_bound = sqrt(2 * omega_debye);
#else
			// use theta(omega - 0.5*|l^2 - k^2|) exactly
			const c_float lower_bound = sqrt(std::max(0.0, bare_dispersion(k) - 2 * omega_debye));
			const c_float upper_bound = sqrt(bare_dispersion(k) + 2 * omega_debye);
#endif

			auto integrand = [this](c_float x) -> c_complex {
				return x * x * sc_expectation_value(x);
				};
			result(u_idx) = Utility::Numerics::Integration::trapezoidal_rule(integrand, lower_bound, upper_bound, N_L);
		}
		constexpr double prefactor = -1. / (2. * M_PI * M_PI);
		result *= static_cast<c_float>(prefactor * U);
		this->Delta.fill_with(result);
		result -= initial_values;
	}

	c_float SCModel::computeCoefficient(const SymbolicOperators::Coefficient& coeff, c_float k, c_float l) const
	{
		if (coeff.name == "\\epsilon_0") 
		{
			if(coeff.momenta.front().momentum_list.front() == 'k')
			{
				return bare_dispersion_to_fermi_level(k);
			}
			else if(coeff.momenta.front().momentum_list.front() == 'l')
			{
				return bare_dispersion_to_fermi_level(l);
			}
			else {
				throw std::invalid_argument("Wrong momentum!");
			}
		} 
		else if (coeff.name == "U")
		{
			if(coeff.momenta[0] == coeff.momenta[1])
			{
				return this->U;
			} 
			else if (coeff.momenta[0].momentum_list.front() == 'k' && coeff.momenta[1].momentum_list.front() == 'l')
			{
				if (omega_debye - abs(bare_dispersion(k) - bare_dispersion(l)) > 0)
				{
					return this->U;
				}
				else 
				{
					return c_float{};
				}
			}
			else {
				throw std::invalid_argument("Wrong momentum!");
			}
		} 
		else 
		{
			throw std::invalid_argument("Coefficient not recognized! " + coeff.name);
		}
	}
	c_complex SCModel::computeTerm(const SymbolicOperators::WickTerm& term, c_float k, c_float l) const
	{
		c_float coeff_value = term.coefficients.empty() ? 1 : this->computeCoefficient(term.coefficients.front(), k, l);
		coeff_value *= term.multiplicity;
	}
}