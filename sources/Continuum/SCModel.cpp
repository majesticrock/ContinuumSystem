#include "SCModel.hpp"
#include <Utility/Numerics/Minimization/Bisection.hpp>
#include <Utility/Numerics/Roots/Bisection.hpp>
#include <algorithm>
#include <numeric>
#include <complex>
//#include <boost/math/tools/roots.hpp>
#include <Utility/Numerics/Roots/Bisection.hpp>
#include <limits>

using Utility::constexprPower;

namespace Continuum {
	SCModel::SCModel(ModelInitializer const& parameters)
		: Delta(2 * DISCRETIZATION, parameters.phonon_coupling* parameters.omega_debye), temperature{ parameters.temperature },
		phonon_coupling{ parameters.phonon_coupling }, omega_debye{ parameters.omega_debye }, fermi_energy{ parameters.fermi_energy },
		fermi_wavevector{ compute_fermiwavevector(fermi_energy) },
		V_OVER_N{ fermi_wavevector > 0 ? 3. * PI * PI / (constexprPower<3>(fermi_wavevector)) : 1 },
		momentumRanges(fermi_wavevector, omega_debye)
	{
		omega_debye += PRECISION;

		auto adjust_kf = [this](c_float kF) {
			fermi_wavevector = kF;
			return bare_dispersion_to_fermi_level(kF) + fock_energy(kF);
		};
		Utility::Numerics::Roots::bisection(adjust_kf, 0.95 * fermi_wavevector, 1.05 * fermi_wavevector, PRECISION, 250);

		Delta = decltype(Delta)::FromAllocator([&](size_t i) -> c_complex {
			const c_float k = momentumRanges.index_to_momentum(i);
			const c_float magnitude = (k < sqrt(2. * (fermi_energy - omega_debye)) || k > sqrt(2. * (fermi_energy + omega_debye)))
				? 0.01 : 0.1;
			if (i < DISCRETIZATION) {
#ifdef _complex
				return magnitude; //std::polar(magnitude, i * 2.0 * PI / (DISCRETIZATION) + 0.5 * PI);
#else
				return std::polar(magnitude, i * 2.0 * PI / (DISCRETIZATION) + 0.5 * PI).real();
#endif
			}
			return (momentumRanges.index_to_momentum(k) > fermi_wavevector ? -0.001 : 0.001 );
			}, 2 * DISCRETIZATION);

		//auto test = [](c_float q){
		//	return 1.0;
		//};
		//std::cout << (I_1(test, 1.)) / PhysicalConstants::em_factor << std::endl;
	}

	std::vector<c_float> SCModel::continuum_boundaries() const
	{
		return {
			2 * this->energy(Utility::Numerics::Minimization::bisection([this](c_float k) { return this->energy(k); },
				momentumRanges.index_to_momentum(0), momentumRanges.index_to_momentum(DISCRETIZATION - 1), PRECISION, 100)),
			2 * this->energy(Utility::Numerics::Minimization::bisection([this](c_float k) { return -this->energy(k); },
				momentumRanges.index_to_momentum(0), momentumRanges.index_to_momentum(DISCRETIZATION - 1), PRECISION, 100))
		};
	}

	c_float SCModel::fock_energy(c_float k) const 
	{
#ifdef _screening
		const c_float ln_factor{ ((1. + _screening * _screening) * fermi_wavevector * fermi_wavevector - k * k) / (4.0 * k * fermi_wavevector) };
		const c_float k_diff{ k - fermi_wavevector };
		const c_float k_sum{ k + fermi_wavevector };

		return -PhysicalConstants::em_factor * fermi_wavevector * 
		(
			1.0 + ln_factor * std::log( (_screening * _screening + k_sum * k_sum) / (_screening * _screening + k_diff * k_diff) )
			+ (_screening / fermi_wavevector) * (std::atan(k_diff / _screening) - std::atan(k_sum / _screening))
		);

#else
		if(is_zero(k - fermi_wavevector)) {
			return -PhysicalConstants::em_factor * fermi_wavevector;
		}

		return -PhysicalConstants::em_factor * fermi_wavevector * (
			1.0 + ((fermi_wavevector * fermi_wavevector - k * k) / (2.0 * k * fermi_wavevector)) 
				* std::log(std::abs((k + fermi_wavevector) / (k - fermi_wavevector)))
		);
#endif
	}

	c_complex SCModel::sc_expectation_value(c_float k) const {
		const auto DELTA = interpolate_delta(k);
		if (is_zero(DELTA)) return 0;
		const c_float E = energy(k);
		if (is_zero(temperature)) {
			return -DELTA / (2 * E);
		}
		return -std::tanh(E / (2 * temperature)) * DELTA / (2 * E);
	}

	c_float SCModel::occupation(c_float k) const {
		const auto DELTA = std::abs(interpolate_delta(k));
		const auto eps_mu = dispersion_to_fermi_level(k);
		if (is_zero(DELTA)) {
			if (is_zero(temperature)) {
				if (is_zero(eps_mu)) return 0.5;
				return (eps_mu < 0 ? 1 : 0);
			}
			return 1. / (1 + std::exp(eps_mu / temperature));
		}
		const c_float E = sqrt(eps_mu * eps_mu + DELTA * DELTA);
		if (is_zero(temperature)) {
			return 0.5 * (1 - (eps_mu / E));
		}
		return 0.5 * (1 - (eps_mu / E) * std::tanh(E / (2 * temperature)));
	}

	void SCModel::iterationStep(const ParameterVector& initial_values, ParameterVector& result) {
		static int step_num = 0;
		result.setZero();
		this->Delta.fill_with(initial_values);
		auto phonon_integrand = [this](c_float x) -> c_complex {
			return x * x * sc_expectation_value(x);
			};
		auto delta_n_wrapper = [this](c_float q) {
				return this->delta_n(q);
				};
		auto sc_wrapper = [this](c_float q) {
			return this->sc_expectation_value(q);
			};

//#pragma omp parallel for
		for (int u_idx = 0; u_idx < DISCRETIZATION; ++u_idx) {
			const c_float k = momentumRanges.index_to_momentum(u_idx);

#ifdef approximate_theta
			// approximate theta(omega - 0.5*|l^2 - k^2|) as theta(omega - eps_k)theta(omega - eps_l)
			if (std::abs(bare_dispersion_to_fermi_level(k)) > omega_debye) {
				continue;
			}
#endif

#ifdef _use_coulomb
#ifdef _screening
			result(u_idx) = integral_screening(sc_wrapper, k);
			result(u_idx + DISCRETIZATION) = integral_screening(delta_n_wrapper, k);
#else
			result(u_idx) = I_1(sc_wrapper, k) + I_2(sc_wrapper, k);
			result(u_idx + DISCRETIZATION) = I_1(delta_n_wrapper, k) + I_2(delta_n_wrapper, k);
#endif
#endif
			//std::cout << g_lower_bound(k) << "\t" << g_upper_bound(k) << std::endl;
			result(u_idx) -= (V_OVER_N * phonon_coupling / (2. * PI * PI))
				* boost::math::quadrature::gauss<double, 30>::integrate( phonon_integrand, g_lower_bound(k), g_upper_bound(k) );
		}

		this->Delta.fill_with(result, 0.5);
		this->Delta.clear_noise(PRECISION);
		result -= initial_values;
		++step_num;
	}

	c_float SCModel::g_lower_bound(c_float k) const
	{
#ifdef approximate_theta
		return sqrt(std::max((static_cast<c_float>(2) * (fermi_energy - omega_debye)), c_float{}));
#else
#ifdef _use_coulomb
		const c_float ALPHA = 2. * omega_debye - phonon_alpha(k);
		auto func = [&](c_float l){
			return this->phonon_beta(l, ALPHA);
			};
		const auto lb = func(momentumRanges.K_MIN);
		const auto ub = func(k);
		if(lb * ub >= c_float{}) return momentumRanges.K_MIN;
		return Utility::Numerics::Roots::bisection(func, momentumRanges.K_MIN, k, 1e-14, 100);
#else
		return std::max(sqrt(k * k - 2. * omega_debye), momentumRanges.K_MIN);
#endif
#endif
	}
	
	c_float SCModel::g_upper_bound(c_float k) const
	{
#ifdef approximate_theta
		return sqrt(2. * (fermi_energy + omega_debye));
#else
#ifdef _use_coulomb
		const c_float ALPHA = 2. * omega_debye + phonon_alpha(k);
		auto func = [&](c_float l){
			return this->phonon_beta(l, -ALPHA);
			};
		const auto lb = func(k);
		const auto ub = func(momentumRanges.K_MAX);
		if(lb * ub >= c_float{}) return momentumRanges.K_MAX;
		return Utility::Numerics::Roots::bisection(func, k, momentumRanges.K_MAX, 1e-14, 100);
#else
		return std::min(sqrt(k * k + 2. * omega_debye), momentumRanges.K_MAX);
#endif
#endif
	}

	c_float SCModel::computeCoefficient(SymbolicOperators::Coefficient const& coeff, c_float first, c_float second) const
	{
		if (coeff.name == "\\epsilon_0")
		{
			return bare_dispersion_to_fermi_level(first);
		}
		else if (coeff.name == "g")
		{
#ifdef approximate_theta
			if (omega_debye > std::abs(bare_dispersion_to_fermi_level(first))
				&& omega_debye > std::abs(bare_dispersion_to_fermi_level(second)))
#else
			if (coeff.momenta[0] == coeff.momenta[1])
			{
				return this->phonon_coupling * this->V_OVER_N;
			}
			if (omega_debye > std::abs(bare_dispersion(first) - bare_dispersion(second)))
#endif
			{
				return this->phonon_coupling * this->V_OVER_N;
			}
			else
			{
				return c_float{};
			}
		} 
		else if(coeff.name == "V") {
#ifdef _screening
			return PhysicalConstants::em_factor / (first * first + _screening * _screening);
#endif
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
			const c_float momentum = momentumRanges.index_to_momentum(k);
			ret.at(SymbolicOperators::Number_Type)[k] = this->occupation(momentum);
			ret.at(SymbolicOperators::SC_Type)[k] = this->sc_expectation_value(momentum);
		}

		return ret;
	}

	c_float SCModel::internal_energy() const 
	{
		auto procedure = [this](c_float k) -> c_float {
			return -k * k * energy(k);
			};
		return (V_OVER_N / (2 * PI * PI)) * boost::math::quadrature::gauss<c_float, 30>::integrate(
			procedure, momentumRanges.K_MIN, momentumRanges.K_MAX);
	}

	c_float SCModel::compute_fermiwavevector(c_float epsilon_F) const
	{
		return PhysicalConstants::em_factor + sqrt((PhysicalConstants::em_factor * PhysicalConstants::em_factor) + 2. * epsilon_F);
	}
}