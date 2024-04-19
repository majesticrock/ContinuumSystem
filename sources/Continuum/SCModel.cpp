#include "SCModel.hpp"
#include "../../../Utility/sources/Numerics/TrapezoidalRule.hpp"

namespace Continuum {
	SCModel::SCModel(ModelInitializer const& parameters)
		: Delta(DISCRETIZATION, parameters.U),  temperature{ parameters.temperature }, U{ parameters.U }, omega_debye{ parameters.omega_debye }
	{ 
	}

	void SCModel::iterationStep(const ParameterVector& initial_values, ParameterVector& result) {
		result.setZero();
		this->Delta.fill_with(initial_values);

		constexpr int N_L = 1000;

		for (int u_idx = 0; u_idx < DISCRETIZATION; ++u_idx) {
			const c_float k = index_to_momentum(u_idx);
			const c_float lower_bound = k*k - omega_debye > 0 ? k*k - omega_debye : 0;
			const c_float upper_bound = k*k + omega_debye;

			auto integrand = [this](c_float k) -> c_complex {
				return k * k * sc_expectation_value(k);
				};
			result(u_idx) = Utility::Numerics::Integration::trapezoidal_rule(integrand, lower_bound, upper_bound, N_L);
		}
		result *= static_cast<c_float>(-4 * M_PI * U);
		this->Delta.fill_with(result);
		result -= initial_values;
	}
}