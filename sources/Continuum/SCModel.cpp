#include "SCModel.hpp"
#define _USE_MATH_DEFINES

namespace Continuum{
        SCModel::SCModel(ModelInitializer const& parameters)
            : temperature{ parameters.temperature }, U{ parameters.U }, omega_debye{ parameters.omega_debye }
        { }

    void SCModel::iterationStep(const ParameterVector& initial_values, ParameterVector& result){
        result.setZero();
        this->Delta.fill_with(initial_values);
        for(int k = 0; k < k + omega_debye_as_index && k < DISCRETIZATION; ++k){
            int l{k - omega_debye_as_index};
            if (l < 0) l = 0;
            for(;l < k + omega_debye_as_index; ++l){
                result(k) += sc_expectation_value(k);
            }
        }
        result *= DELTA_K;

        result *= static_cast<c_float>(-4 * M_PI * U);
        this->Delta.fill_with(result);
        result -= initial_values;
    }
}