#include "SplineContainer.hpp"

namespace Continuum {
	SplineContainer::SplineContainer(const std::vector<c_complex>& ys, c_float begin,
			c_float lower_step, c_float middle_step, c_float upper_step, int n_lower, int n_middle)
		: construct{begin, lower_step, middle_step, upper_step, 
			begin + lower_step * n_lower, begin + lower_step * n_lower + middle_step * n_middle, 
			n_lower, n_middle, static_cast<int>(ys.size()) - n_lower - n_middle},
		lower_spline(ys.data(), n_lower + 1, begin, lower_step),
		middle_spline(ys.data() + n_lower, n_middle + 1, construct.end_lower, middle_step),
		upper_spline(ys.data() + n_lower + n_middle, construct.n_upper, construct.end_middle, upper_step)
	{}

	void SplineContainer::set_new_ys(const std::vector<c_complex>& ys)
	{
		lower_spline = Spline(ys.data(), construct.n_lower + 1, construct.begin, construct.lower_step);
		middle_spline = Spline(ys.data() + construct.n_lower, construct.n_middle + 1, construct.end_lower, construct.middle_step);
		upper_spline = Spline(ys.data() + construct.n_lower + construct.n_middle, construct.n_upper, construct.end_middle, construct.upper_step);
	}
}