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
#include <mrock/utility/ComplexNumberIterators.hpp>

namespace Continuum {
	template <class T>
	struct ComplexSpline {
		using Spline = boost::math::interpolators::_SPLINE_TYPE <T>;
		Spline real, imag;
		ComplexSpline(std::complex<T> const* const data, int length, T left_endpoint, T step_size)
		{
			auto real_begin = mrock::utility::make_real_part_iterator(data);
			auto imag_begin = mrock::utility::make_imag_part_iterator(data);
			auto real_end = mrock::utility::make_real_part_iterator_end(data, length);
			auto imag_end = mrock::utility::make_imag_part_iterator_end(data, length);
			real = Spline(real_begin, real_end, left_endpoint, step_size);
			imag = Spline(imag_begin, imag_end, left_endpoint, step_size);
		}
		ComplexSpline() = default;

		inline std::complex<T> operator()(T k) const {
			return std::complex<T>(real(k), imag(k));
		}
	};

	class SplineContainer {
	private:
		// used for construction of the splines with new ys (and same xs)
		struct _construct {
			c_float begin;
			c_float lower_step, middle_step, upper_step;
			c_float end_lower, end_middle;
			int n_lower, n_middle, n_upper;
		};
		_construct construct;

		using Spline =
#ifdef _complex
			ComplexSpline<c_float>;
#else
			boost::math::interpolators::_SPLINE_TYPE <c_float>;
#endif
		Spline lower_spline, middle_spline, upper_spline;

	public:
		SplineContainer(const std::vector<c_complex>& ys, c_float begin,
			c_float lower_step, c_float middle_step, c_float upper_step,
			int n_lower, int n_middle);
		SplineContainer() = default;

		void set_new_ys(const std::vector<c_complex>& ys);

		inline c_complex operator()(c_float k) const {
			if (k < construct.end_lower) {
				return lower_spline(k);
			}
			if (k < construct.end_middle) {
				return middle_spline(k);
			}
			return upper_spline(k);
		}
	};
}

#undef _quintic_spline
#undef _SPLINE_TYPE