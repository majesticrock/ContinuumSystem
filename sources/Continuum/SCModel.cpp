#include "SCModel.hpp"
#include "../../../Utility/sources/Numerics/Integration/TrapezoidalRule.hpp"
#include "../../../Utility/sources/ConstexprPower.hpp"
#include "../../../Utility/sources/Numerics/Minimization/Bisection.hpp"
#include <algorithm>
#include <numeric>

using Utility::constexprPower;

namespace Continuum {
	SCModel::SCModel(ModelInitializer const& parameters)
		: Delta(DISCRETIZATION, parameters.U* parameters.omega_debye), temperature{ parameters.temperature },
		U{ parameters.U }, omega_debye{ parameters.omega_debye }, fermi_energy{ parameters.fermi_energy },
		fermi_wavevector{ sqrt(2 * parameters.fermi_energy) },
		V_OVER_N{ fermi_wavevector > 0 ? 3. * M_PI * M_PI / (constexprPower<3>(fermi_wavevector)) : 1 },
		U_MAX{ sqrt(2 * (fermi_energy + omega_debye)) - fermi_wavevector },
		U_MIN{ fermi_energy > omega_debye ? sqrt(2 * (fermi_energy - omega_debye)) - fermi_wavevector : -fermi_wavevector },
		STEP{ (U_MAX - U_MIN) / (DISCRETIZATION - 1) }
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
		const auto DELTA = std::norm(interpolate_delta(k));
		const auto eps_mu = bare_dispersion_to_fermi_level(k);
		if (is_zero(DELTA)) {
			if (is_zero(temperature)) {
				if (is_zero(eps_mu)) return 0.5;
				return (eps_mu < 0 ? 1 : 0);
			}
			return 1. / (1 + std::exp(eps_mu / temperature));
		}
		const c_float E = sqrt(eps_mu * eps_mu + DELTA);
		if (is_zero(temperature)) {
			return 0.5 * (1 - (eps_mu / E));
		}
		return 0.5 * (1 - (eps_mu / E) * std::tanh(E / (2 * temperature)));
	}

	void SCModel::iterationStep(const ParameterVector& initial_values, ParameterVector& result) {
		result.setZero();
		this->Delta.fill_with(initial_values);

		constexpr int N_L = 1000;

		for (int u_idx = 0; u_idx < DISCRETIZATION; ++u_idx) {
			const c_float k = index_to_momentum(u_idx);
#ifdef approximate_theta
			// approximate theta(omega - 0.5*|l^2 - k^2|) as theta(omega - eps_k)theta(omega - eps_l)
			if (std::abs(bare_dispersion_to_fermi_level(k)) > omega_debye) {
				continue;
			}
#endif
			auto integrand = [this](c_float x) -> c_complex {
				return x * x * sc_expectation_value(x);
				};
			result(u_idx) = Utility::Numerics::Integration::trapezoidal_rule(integrand, u_lower_bound(k), u_upper_bound(k), N_L);
		}
		const c_float prefactor = -V_OVER_N / (2. * M_PI * M_PI);
		result *= prefactor * U;
		this->Delta.fill_with(result);
		result -= initial_values;
	}

	c_float SCModel::computeCoefficient(SymbolicOperators::Coefficient const& coeff, c_float first, c_float second /*=c_float{}*/) const
	{
		if (coeff.name == "\\epsilon_0")
		{
			return bare_dispersion_to_fermi_level(first);
		}
		else if (coeff.name == "U")
		{
			if (coeff.momenta[0] == coeff.momenta[1])
			{
				return this->U * this->V_OVER_N;
			}
#ifdef approximate_theta
			if (omega_debye > std::abs(bare_dispersion_to_fermi_level(first))
				&& omega_debye > std::abs(bare_dispersion_to_fermi_level(second)))
#else
			if (omega_debye > std::abs(bare_dispersion(first) - bare_dispersion(second)))
#endif
			{
				return this->U * this->V_OVER_N;
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