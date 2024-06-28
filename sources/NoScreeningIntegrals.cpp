		static constexpr c_float OFF_SET_FACTOR = 0.03;
		template<class ExpectationValues>
		inline auto I_1(ExpectationValues const& expecs, c_float k) const {
			auto integrand = [&expecs, &k](c_float x) {
				const c_float tanh_x = std::tanh(0.5 * x);
				return expecs(k * tanh_x) * x * tanh_x / Utility::constexprPower<2>(std::cosh(0.5 * x));
				};

			const c_float alpha{ momentumRanges.K_MIN / k };//
			if(is_zero(1. - alpha)) return decltype(expecs(k)){};
			const c_float lower_bound = std::log( (1. + alpha) / (1. - alpha) );
			constexpr c_float upper_bound = 30; // a priori error estimation of machine epsilon

			const c_float prefactor = 0.5 * PhysicalConstants::em_factor * k;
			constexpr c_float CUT_OFF_CONSTANT = 1.16034524813597e-11;//1.73136903588109e-7;

			if(k > fermi_wavevector){
				const c_float middle_bound = 2. * std::atanh(fermi_wavevector / k);
				const c_float k_f_offset = OFF_SET_FACTOR * DISCRETIZATION  * momentumRanges.STEP;

				if(k > fermi_wavevector + k_f_offset){
					const c_float middle_middle_bound = 2. * std::atanh((fermi_wavevector + k_f_offset) / k);

					return prefactor * (CUT_OFF_CONSTANT * expecs(k)
						+ boost::math::quadrature::gauss<double, n_gauss>::integrate(integrand, lower_bound, middle_bound)
						+ boost::math::quadrature::gauss<double, n_gauss>::integrate(integrand, middle_bound, middle_middle_bound)
						+ boost::math::quadrature::gauss<double, n_gauss>::integrate(integrand, middle_middle_bound, upper_bound));
				}

				return prefactor * (CUT_OFF_CONSTANT * expecs(k)
					+ boost::math::quadrature::gauss<double, n_gauss>::integrate(integrand, lower_bound, middle_bound)
					+ boost::math::quadrature::gauss<double, n_gauss>::integrate(integrand, middle_bound, upper_bound));
			}

			return prefactor * (CUT_OFF_CONSTANT * expecs(k)
				+ boost::math::quadrature::gauss<double, n_gauss>::integrate(integrand, lower_bound, upper_bound));
		}

		template<class ExpectationValues>
		inline auto I_2(ExpectationValues const& expecs, c_float k) const {
			auto integrand = [&expecs, &k](c_float x) {
				const c_float coth_x = 1. / std::tanh(0.5 * x);
				return expecs(k * coth_x) * x * coth_x / Utility::constexprPower<2>(std::sinh(0.5 * x));
				};
			const c_float beta{ momentumRanges.K_MAX / k };
			if(is_zero(beta - 1.)) return decltype(expecs(k)){};
			const c_float lower_bound = std::log( (beta + 1.) / (beta - 1.) );
			constexpr c_float upper_bound = 30; // a priori error estimation of machine epsilon

			const c_float prefactor = 0.5 * PhysicalConstants::em_factor * k;
			constexpr c_float CUT_OFF_CONSTANT = 1.16034524813640e-11;//1.73136904981569e-7;

			if(k < fermi_wavevector){
				const c_float middle_bound = 2. * std::atanh(k / fermi_wavevector); // acoth(x) = atanh(1/x)
				const c_float k_f_offset = OFF_SET_FACTOR * DISCRETIZATION * momentumRanges.STEP;

				if(k < fermi_wavevector - k_f_offset){
					const c_float middle_middle_bound = 2. * std::atanh(k / (fermi_wavevector - k_f_offset));

					return prefactor * (CUT_OFF_CONSTANT * expecs(k)
						+ boost::math::quadrature::gauss<double, n_gauss>::integrate(integrand, lower_bound, middle_bound)
						+ boost::math::quadrature::gauss<double, n_gauss>::integrate(integrand, middle_bound, middle_middle_bound)
						+ boost::math::quadrature::gauss<double, n_gauss>::integrate(integrand, middle_middle_bound, upper_bound));
				}

				return prefactor * (CUT_OFF_CONSTANT * expecs(k)
					+ boost::math::quadrature::gauss<double, n_gauss>::integrate(integrand, lower_bound, middle_bound)
					+ boost::math::quadrature::gauss<double, n_gauss>::integrate(integrand, middle_bound, upper_bound));
			}

			return prefactor * (CUT_OFF_CONSTANT * expecs(k)
				+ boost::math::quadrature::gauss<double, n_gauss>::integrate(integrand, lower_bound, upper_bound));
		}