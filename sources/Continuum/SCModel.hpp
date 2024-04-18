#pragma once
#include "GlobalDefinitions.hpp"
#include "ModelAttributes.hpp"
#include <cmath>
#include <limits>
#include <boost/math/special_functions/pow.hpp>

namespace Continuum{
    struct ModelInitializer {
        c_float temperature;
        c_float U;
        c_float omega_debye;
    };

    class SCModel {
    public:
        ModelAttributes<c_complex> Delta;

        c_float temperature{};
        c_float U{};
        c_float omega_debye{};

        static constexpr c_float CUT_OFF = std::numeric_limits<c_float>::epsilon();
        
        inline c_float energy(int k) const {
            return sqrt(boost::math::pow<4>(index_to_momentum(k)) + Delta[k].norm());
        }
        inline c_float sc_expectation_value(int k) const {
            const c_float E = energy(k);
            if( is_zero(E) ) return 0;
            if( is_zero(temperature) ){
                return Delta[k] / (2 * E);
            }
            return std::tanh(E / (2 * temperature)) * Delta[k] / (2 * E);
            
        };
        void iterationStep(const ParameterVector& initial_values, ParamterVector& result);

        SCModel(ModelInitializer const& parameters);
    };
}