#pragma once
//#define _quintic_spline
#ifdef _quintic_spline // quintic does not work because it has no default constructor, which we need
#include <boost/math/interpolators/cardinal_quintic_b_spline.hpp>
#define _SPLINE_TYPE cardinal_quintic_b_spline
#else
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
#define _SPLINE_TYPE cardinal_cubic_b_spline
#endif

#include "GlobalDefinitions.hpp"
#include <vector>

namespace Continuum {
    class SplineContainer {
    private:
    // used for construction of the splines with new ys (and same xs)
        struct _construct{
            c_float begin;
            c_float lower_step, middle_step, upper_step;
            c_float end_lower, end_middle;
            int n_lower, n_middle, n_upper;
        };
        _construct construct;
        
        using Spline = boost::math::interpolators:: _SPLINE_TYPE <c_float>;
        Spline lower_spline, middle_spline, upper_spline;

    public:
        SplineContainer(const std::vector<double>& ys, c_float begin, 
            c_float lower_step, c_float middle_step, c_float upper_step, 
            int n_lower, int n_middle);
        SplineContainer() = default;

        void set_new_ys(const std::vector<double>& ys);

        inline c_float operator()(c_float k) const {
            if(k < construct.end_lower) {
                return lower_spline(k);
            }
            if(k < construct.end_middle) {
                return middle_spline(k);
            }
            return upper_spline(k);
        }
    };
}

#undef _quintic_spline
#undef _SPLINE_TYPE