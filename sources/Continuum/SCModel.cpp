#include "SCModel.hpp"
#include <Utility/Numerics/Minimization/Bisection.hpp>
#include <Utility/Selfconsistency/BroydenSolver.hpp>
#include <algorithm>
#include <numeric>
#include <complex>
#include <Utility/Numerics/Roots/Bisection.hpp>
#include <Utility/better_to_string.hpp>
#include <limits>

using Utility::constexprPower;

namespace Continuum {
	SCModel::SCModel(ModelInitializer const& parameters)
		: Delta(2 * DISCRETIZATION + 1, c_complex{}), 
		temperature{ parameters.temperature }, phonon_coupling{ parameters.phonon_coupling }, 
		omega_debye{ parameters.omega_debye }, fermi_energy{ parameters.fermi_energy },
		 screening{ parameters.screening },coulomb_scaling{ parameters.coulomb_scaling },
		fermi_wavevector{ parameters.fermi_wavevector },
		momentumRanges(&fermi_wavevector, omega_debye)
	{
		std::cout << "Fock(k_F) = " << fock_energy(fermi_wavevector) << "  xi(k_F) = " << dispersion_to_fermi_level(fermi_wavevector) << std::endl;
		Delta = decltype(Delta)::FromAllocator([&](size_t i) -> c_complex {
			const c_float k = momentumRanges.index_to_momentum(i);
			const c_float magnitude = (k < sqrt(2. * (fermi_energy - omega_debye)) || k > sqrt(2. * (fermi_energy + omega_debye))) ? 0.01 : 0.1;
			if (i < DISCRETIZATION) {
#ifdef _complex
				return std::polar(magnitude, PI * (static_cast<double>(i) / static_cast<double>(DISCRETIZATION)));
#else
				return magnitude * std::cos(PI * (static_cast<double>(i) / static_cast<double>(DISCRETIZATION) - 0.5));
#endif
			}
			if (i == 2 * DISCRETIZATION) return c_complex{};
			return c_complex{};
			}, 2 * DISCRETIZATION + 1);
		set_splines();
	}

	void SCModel::set_new_parameters(ModelInitializer const& parameters) 
	{
		this->temperature = parameters.temperature; 
		this->omega_debye = parameters.omega_debye; 
		this->fermi_energy = parameters.fermi_energy; 
		this->phonon_coupling = parameters.phonon_coupling; 
		this->coulomb_scaling = parameters.coulomb_scaling; 
		this->fermi_wavevector = parameters.fermi_wavevector;
		try {
			this->momentumRanges = MomentumRanges(&this->fermi_wavevector, parameters.omega_debye);
		} 
		catch(...) {
			std::cerr << parameters << std::endl;
			throw;
		}
		if(is_zero(Delta[DISCRETIZATION / 2])) {
			Delta = decltype(Delta)::FromAllocator([&](size_t i) -> c_complex {
			const c_float k = momentumRanges.index_to_momentum(i);
			const c_float magnitude = (k < sqrt(2. * (fermi_energy - omega_debye)) || k > sqrt(2. * (fermi_energy + omega_debye))) ? 0.01 : 0.1;
			if (i < DISCRETIZATION) {
#ifdef _complex
				return std::polar(magnitude, 1.3 + PI * (static_cast<double>(i) / static_cast<double>(DISCRETIZATION)));
#else
				return magnitude * std::cos(PI * (static_cast<double>(i) / static_cast<double>(DISCRETIZATION) - 0.5));
#endif
			}
			if (i == 2 * DISCRETIZATION) return c_complex{};
			return (momentumRanges.index_to_momentum(k) > fermi_wavevector ? -0.001 : 0.001 );
			}, 2 * DISCRETIZATION + 1);
		}

		set_splines();
		auto solver = Utility::Selfconsistency::make_broyden<c_complex>(this, &Delta, 200);
		solver.compute(true);
	}

	c_complex SCModel::k_infinity_integral() const
	{
		auto integrand = [this](c_float q) {
			return q * q * this->sc_expectation_value(q);
		};
		const c_float prefactor = 1. + 2. * PhysicalConstants::em_factor * coulomb_scaling / momentumRanges.K_MAX;
		return momentumRanges.integrate(integrand) / prefactor;
	}

	c_complex SCModel::k_zero_integral() const
	{
		auto integrand = [this](c_float q) {
			return (q * q / (screening * screening + q * q)) * this->sc_expectation_value(q);
		};
		const c_float prefactor = 2. * coulomb_scaling * PhysicalConstants::em_factor;
		return momentumRanges.integrate(integrand) * prefactor;
	}

	c_float SCModel::fock_energy(c_float k) const 
	{
#ifdef screening
		if(is_zero(k)) {
			return -coulomb_scaling * PhysicalConstants::em_factor * fermi_wavevector * (
				3.0 - 2.0 * (screening / fermi_wavevector) * std::atan(fermi_wavevector / screening)
			);
		}

		const c_float k_diff{ k - fermi_wavevector };
		const c_float k_sum{ k + fermi_wavevector };
		const c_float ln_factor{ (screening * screening + fermi_wavevector * fermi_wavevector - k * k) / (2.0 * k * fermi_wavevector) };
		return -coulomb_scaling * PhysicalConstants::em_factor * fermi_wavevector * 
		(
			1.0 + ln_factor * log_expression(k_sum, k_diff) 
			+ (screening / fermi_wavevector) * (std::atan(k_diff / screening) - std::atan(k_sum / screening))
		);

#else
		if(is_zero(k - fermi_wavevector)) {
			return -coulomb_scaling * PhysicalConstants::em_factor * fermi_wavevector;
		}

		return -coulomb_scaling * PhysicalConstants::em_factor * fermi_wavevector * (
			1.0 + ((fermi_wavevector * fermi_wavevector - k * k) / (2.0 * k * fermi_wavevector)) 
				* std::log(std::abs((k + fermi_wavevector) / (k - fermi_wavevector)))
		);
#endif
	}

	c_complex SCModel::sc_expectation_value_index(int k) const
	{
		if (is_zero(Delta[k])) return c_complex{};
		const c_float E = energy_index(k);
		if (is_zero(temperature)) {
			return -Delta[k] / (2 * E);
		}
		return -std::tanh(E / (2 * temperature)) * Delta[k] / (2 * E);
	}

	c_float SCModel::occupation_index(int k) const
	{
		const auto eps_mu = dispersion_to_fermi_level_index(k);
		if (is_zero(Delta[k])) {
			if (is_zero(temperature)) {
				if (is_zero(eps_mu)) return 0.5;
				return (eps_mu < 0 ? 1 : 0);
			}
			return 1. / (1 + std::exp(eps_mu / temperature));
		}
		const c_float E = sqrt(eps_mu * eps_mu + std::norm(Delta[k]));
		if (is_zero(temperature)) {
			return 0.5 * (1 - (eps_mu / E));
		}
		return 0.5 * (1 - (eps_mu / E) * std::tanh(E / (2 * temperature)));
	}

	void SCModel::iterationStep(const ParameterVector& initial_values, ParameterVector& result) {
		static int step_num = 0;
		result.setZero();
		this->Delta.fill_with(initial_values);
		this->get_expectation_values();
		this->occupation.set_new_ys(_expecs[SymbolicOperators::Number_Type]);
		this->sc_expectation_value.set_new_ys(_expecs[SymbolicOperators::SC_Type]);

		auto delta_n_wrapper = [this](c_float q) {
				return this->delta_n(q);
				};

//#pragma omp parallel for
		for (MomentumIterator it(&momentumRanges); it < DISCRETIZATION; ++it) {
			result(it.idx) = integral_screening(sc_expectation_value, it.k);
#ifndef mielke_coulomb
			result(it.idx + DISCRETIZATION) = integral_screening(delta_n_wrapper, it.k);
#endif
#ifdef approximate_theta
			// approximate theta(omega - 0.5*|l^2 - k^2|) as theta(omega - eps_k)theta(omega - eps_l)
			if (std::abs(phonon_alpha(it.k) - 2. * fermi_energy) > 2.0 * omega_debye) {
				continue;
			}
#endif
			result(it.idx) -= integral_phonon(sc_expectation_value, it.k);
		}
		result(2 * DISCRETIZATION) = k_infinity_integral();

		this->Delta.fill_with(result, 0.5);
		this->Delta.clear_noise(PRECISION);
		result -= initial_values;
		++step_num;
	}

	c_float SCModel::g_lower_bound(c_float k) const
	{
#ifdef approximate_theta
		const c_float ALPHA = 2. * fermi_energy - 2. * omega_debye;
#else
		const c_float ALPHA = phonon_alpha(k) - 2. * omega_debye;
#endif
		auto func = [&](c_float l){
			return this->phonon_beta(l, ALPHA);
			};
#ifdef approximate_theta
		return Utility::Numerics::Roots::bisection(func, momentumRanges.K_MIN, fermi_wavevector, PRECISION, 200);
#else
		const auto lb = func(momentumRanges.K_MIN);
		const auto ub = func(k);
		if(lb * ub >= c_float{}) return momentumRanges.K_MIN;
		return Utility::Numerics::Roots::bisection(func, momentumRanges.K_MIN, k, PRECISION, 200);
#endif
	}
	
	c_float SCModel::g_upper_bound(c_float k) const
	{
#ifdef approximate_theta
		const c_float ALPHA = 2. * fermi_energy + 2. * omega_debye;
#else
		const c_float ALPHA = phonon_alpha(k) + 2. * omega_debye;
#endif
		auto func = [&](c_float l){
			return this->phonon_beta(l, ALPHA);
			};
#ifdef approximate_theta
		return Utility::Numerics::Roots::bisection(func, fermi_wavevector, momentumRanges.K_MAX, PRECISION, 200);
#else
		const auto lb = func(k);
		const auto ub = func(momentumRanges.K_MAX);
		if(lb * ub >= c_float{}) return momentumRanges.K_MAX;
		return Utility::Numerics::Roots::bisection(func, k, momentumRanges.K_MAX, PRECISION, 200);
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
			if (omega_debye > std::abs(bare_dispersion_to_fermi_level(first) + fock_energy(first))
				&& omega_debye > std::abs(bare_dispersion_to_fermi_level(second) + fock_energy(second)))
#else
			if (2. * omega_debye > std::abs(phonon_alpha(first) - phonon_alpha(second)))
#endif
			{
				return this->phonon_coupling;
			}
			else
			{
				return c_float{};
			}
		} 
		else if(coeff.name == "V") {
			if(coeff.momenta.front().is_zero())
				return (coulomb_scaling / PhysicalConstants::vacuum_permitivity) / (screening * screening);
			return coulomb_scaling / (2 * first * second * PhysicalConstants::vacuum_permitivity) * log_expression(first + second, first - second);
		}
		else
		{
			throw std::invalid_argument("Coefficient not recognized! " + coeff.name);
		}
	}

	c_float SCModel::phonon_alpha(const c_float k) const {
		return 2. * bare_dispersion(k) + 2. * fock_energy(k);
	}

	const std::map<SymbolicOperators::OperatorType, std::vector<c_complex>>& SCModel::get_expectation_values() const
	{
		if(_expecs.empty()) {
			_expecs.emplace(SymbolicOperators::Number_Type, std::vector<c_complex>(DISCRETIZATION));
			_expecs.emplace(SymbolicOperators::SC_Type, std::vector<c_complex>(DISCRETIZATION));
		}
		for (int k = 0; k < DISCRETIZATION; ++k) {
			_expecs.at(SymbolicOperators::Number_Type)[k] = this->occupation_index(k);
			_expecs.at(SymbolicOperators::SC_Type)[k] = this->sc_expectation_value_index(k);
		}

		return _expecs;
	}

	c_complex SCModel::interpolate_delta(c_float k) const {
		const int index = momentumRanges.momentum_to_floor_index(k);
		if (index >= DISCRETIZATION - 1) // Assuming Delta(k) = 0 for k -> infinity
			return (index >= DISCRETIZATION ? c_complex{} : Delta[DISCRETIZATION - 1]);
		if (index < 0) // Assuming Delta(k) = 0 for k->0
			return c_complex{};
		return Utility::Numerics::interpolate_from_vector<n_interpolate>(k, momentumRanges, Delta, shifted_index(index));
	}

	c_float SCModel::interpolate_delta_n(c_float k) const {
		const int index = momentumRanges.momentum_to_floor_index(k);
		if (index >= DISCRETIZATION - 1) // Assuming Delta(k) = 0 for k -> infinity
			return (index >= DISCRETIZATION ? c_float{} : std::real(Delta[2 * DISCRETIZATION - 1]));
		if (index < 0) // Assuming Delta(k) = const for k->0
			return c_float{};
		return std::real(Utility::Numerics::interpolate_from_vector<n_interpolate>(k, momentumRanges, Delta, shifted_index(index), DISCRETIZATION));
	}

	c_float SCModel::internal_energy() const 
	{
		if(! is_zero(Delta[DISCRETIZATION / 2])) {
			// Order: -62.0978727910253 eV
			auto procedure = [this](c_float k) -> c_float {
				return -k * k * energy(k) / constexprPower<3>(momentumRanges.K_MAX);
				};
			return momentumRanges.integrate(procedure);
		}
		// No order: -2.25329965038867 eV * 2 because of the spin
		auto procedure = [this](c_float k) -> c_float {
			return k * k * dispersion_to_fermi_level(k) / constexprPower<3>(momentumRanges.K_MAX);
			};
		return 2. * momentumRanges.integrate(procedure, momentumRanges.K_MIN, fermi_wavevector);
	}

	std::string SCModel::info() const {
		return "SCModel // T=" + std::to_string(temperature) + "K   g="
			+ std::to_string(phonon_coupling) + "eV   omega_D=" 
			+ std::to_string(1e3 * omega_debye) + "meV   E_F="
			+ std::to_string(fermi_energy) + "eV   alpha=" + std::to_string(coulomb_scaling);
			+ "lambda=" + std::to_string(screening) + "sqrt(eV)";
	}

	std::string SCModel::to_folder() const {
		auto improved_string = [](c_float number) -> std::string {
			if (std::floor(number) == number) {
				// If the number is a whole number, format it with one decimal place
				std::ostringstream out;
				out.precision(1);
				out << std::fixed << number;
				return out.str();
			}
			else {
				std::string str = Utility::better_to_string(number, std::chars_format::fixed);
				// Remove trailing zeroes
				str.erase(str.find_last_not_of('0') + 1, std::string::npos);
				str.erase(str.find_last_not_of('.') + 1, std::string::npos);
				return str;
			}
			};

		return "T=" + improved_string(temperature) 
			+ "/coulomb_scaling=" + improved_string(coulomb_scaling)
			+ "/screening=" + improved_string(screening)
			+ "/k_F=" + improved_string(fermi_wavevector)
			+ "/g=" + improved_string(phonon_coupling) 
			+ "/omega_D=" + improved_string(1e3 * omega_debye) + "/";		
	}

	std::vector<c_complex> SCModel::phonon_gap() const
	{
		std::vector<c_complex> ret(DISCRETIZATION);
		for(MomentumIterator it(&momentumRanges); it < DISCRETIZATION; ++it) {
#ifdef approximate_theta
			if (std::abs(phonon_alpha(it.k) - 2. * fermi_energy) > 2.0 * omega_debye) {
				continue;
			}
#endif
			ret[it.idx] = -integral_phonon(sc_expectation_value, it.k);
		}
		return ret;
	}

	std::vector<c_complex> SCModel::coulomb_gap() const
	{
		std::vector<c_complex> ret(DISCRETIZATION);
		for(MomentumIterator it(&momentumRanges); it < DISCRETIZATION; ++it) {
			ret[it.idx] = integral_screening([this](c_float q) { return this->sc_expectation_value(q); }, it.k );
		}
		return ret;
	}

	std::vector<c_float> SCModel::single_particle_dispersion() const 
	{
		std::vector<c_float> ret(DISCRETIZATION);
		for(MomentumIterator it(&momentumRanges); it < DISCRETIZATION; ++it) {
			ret[it.idx] = bare_dispersion_to_fermi_level(it.k) + fock_energy(it.k);
		}
		return ret;
	}

	void SCModel::set_splines() 
	{
		this->get_expectation_values();
		this->occupation = SplineContainer(_expecs[SymbolicOperators::Number_Type], momentumRanges.K_MIN,
			momentumRanges.LOWER_STEP, momentumRanges.INNER_STEP, momentumRanges.UPPER_STEP,
			_OUTER_DISC, _INNER_DISC);
		this->sc_expectation_value = SplineContainer(_expecs[SymbolicOperators::SC_Type], momentumRanges.K_MIN,
			momentumRanges.LOWER_STEP, momentumRanges.INNER_STEP, momentumRanges.UPPER_STEP,
			_OUTER_DISC, _INNER_DISC);
	}
}