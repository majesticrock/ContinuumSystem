#include "SCModel.hpp"

namespace Continuum {
	SCModel::SCModel(ModelInitializer const& parameters)
		: temperature{ parameters.temperature }, U{ parameters.U }, omega_debye{ parameters.omega_debye }
	{ }

	void SCModel::iterationStep(const ParameterVector& initial_values, ParameterVector& result) {
		result.setZero();
		this->Delta.fill_with(initial_values);

		constexpr int N_L = 1000;

		for (int u_idx = 0; u_idx < DISCRETIZATION; ++u_idx) {
			const c_float u = STEP * u_idx;
			const c_float k = index_to_momentum(u_idx);
			const c_float lower_bound = k - omega_debye > 0 ? k - omega_debye : 0;
			const c_float upper_bound = k + omega_debye;
			const c_float integral_range = upper_bound - lower_bound;
			const c_float l_step = integral_range / N_L;

			result(u_idx) = 0.5 * (lower_bound * lower_bound * sc_expectation_value(lower_bound)
				+ upper_bound * upper_bound * sc_expectation_value(upper_bound));
			for (int l_idx = 1; l_idx < N_L; ++l_idx) {
				const c_float l = lower_bound + l_idx * l_step;
				result(u_idx) += l * l * sc_expectation_value(l);
			}
			result(u_idx) *= STEP;
		}

		result *= STEP;

		result *= static_cast<c_float>(-4 * M_PI * U);
		this->Delta.fill_with(result);
		result -= initial_values;
	}
}