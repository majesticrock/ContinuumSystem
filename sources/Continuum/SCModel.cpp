#include "SCModel.hpp"
#include "../../../Utility/sources/Numerics/TrapezoidalRule.hpp"
#include <algorithm>

//#define approximate_theta

namespace Continuum {
	SCModel::SCModel(ModelInitializer const& parameters)
		: Delta(DISCRETIZATION, parameters.U),  temperature{ parameters.temperature }, U{ parameters.U }, 
		omega_debye{ parameters.omega_debye }, chemical_potential{ parameters.chemical_potential }
	{ 
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
}