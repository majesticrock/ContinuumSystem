#include "SCModel.hpp"
#include <mrock/utility/Numerics/Minimization/Bisection.hpp>
#include <mrock/utility/Selfconsistency/BroydenSolver.hpp>
#include <algorithm>
#include <numeric>
#include <complex>
#include <mrock/utility/Numerics/Roots/Bisection.hpp>
#include <mrock/utility/better_to_string.hpp>
#include <limits>

using mrock::utility::constexprPower;

namespace Continuum {
	SCModel::SCModel(ModelInitializer const& parameters)
		: Delta(2 * DISCRETIZATION + 1, c_complex{}),
		temperature{ parameters.temperature }, phonon_coupling{ parameters.phonon_coupling },
		omega_debye{ parameters.omega_debye }, fermi_energy{ parameters.fermi_energy },
		screening_ratio{ parameters.screening_ratio }, screening{ parameters.screening },
		coulomb_scaling{ parameters.coulomb_scaling },
		fermi_wavevector{ parameters.fermi_wavevector }, rho_F{ parameters.rho_F },
		momentumRanges(&fermi_wavevector, omega_debye, parameters.x_cut)
	{
		phonon_interaction.set_parent(this);
		std::cout << "Fock(k_F) = " << __fock_coulomb(fermi_wavevector) << "  xi(k_F) = " << dispersion_to_fermi_level(fermi_wavevector) << std::endl;
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
		this->screening_ratio = parameters.screening_ratio;
		this->screening = parameters.screening;
		this->fermi_wavevector = parameters.fermi_wavevector;
		this->rho_F = parameters.rho_F;
		try {
			this->momentumRanges = MomentumRanges(&this->fermi_wavevector, parameters.omega_debye, parameters.x_cut);
		}
		catch (...) {
			std::cerr << parameters << std::endl;
			throw;
		}
		if (is_zero(Delta[DISCRETIZATION / 2])) {
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
				return (momentumRanges.index_to_momentum(k) > fermi_wavevector ? -0.001 : 0.001);
				}, 2 * DISCRETIZATION + 1);
		}
		std::cout << "Working on " << info() << std::endl;
		set_splines();
#ifdef _iterative_selfconsistency
		auto solver = mrock::utility::Selfconsistency::make_iterative<c_complex>(this, &Delta);
#else
		auto solver = mrock::utility::Selfconsistency::make_broyden<c_complex>(this, &Delta, 200);
#endif
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

	c_float SCModel::fock_coulomb(c_float k) const
	{
		if (is_zero(coulomb_scaling)) return c_float{};
		if (is_zero(k)) {
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
	}

	c_complex SCModel::sc_expectation_value_index(int k) const
	{
		if (is_zero(Delta[k])) return c_complex{};
		const c_float E = energy_index(k);
		if (is_zero(temperature)) {
			return -Delta[k] / (2 * E);
		}
		return -std::tanh(E / (2. * PhysicalConstants::k_B * temperature)) * Delta[k] / (2 * E);
	}

	c_float SCModel::occupation_index(int k) const
	{
		const auto eps_mu = dispersion_to_fermi_level_index(k);
		if (is_zero(Delta[k])) {
			if (is_zero(temperature)) {
				if (is_zero(eps_mu)) return 0.5;
				return (eps_mu < 0 ? 1 : 0);
			}
			return 1. / (1 + std::exp(eps_mu / (PhysicalConstants::k_B * temperature)));
		}
		const c_float E = sqrt(eps_mu * eps_mu + std::norm(Delta[k]));
		if (is_zero(temperature)) {
			return 0.5 * (1 - (eps_mu / E));
		}
		return 0.5 * (1 - (eps_mu / E) * std::tanh(E / (2. * PhysicalConstants::k_B * temperature)));
	}

	void SCModel::iteration_step(const ParameterVector& initial_values, ParameterVector& result) {
		static int step_num = 0;
		result.setZero();
		this->Delta.fill_with(initial_values);
		this->get_expectation_values();
		this->occupation.set_new_ys(_expecs[mrock::symbolic_operators::Number_Type]);
		this->sc_expectation_value.set_new_ys(_expecs[mrock::symbolic_operators::SC_Type]);

#if !defined(NO_FOCK_COULOMB) || !defined(NO_FOCK_PHONON) 
		auto delta_n_wrapper = [this](c_float q) {
			return this->delta_n(q);
			};
#endif
		//#pragma omp parallel for
		for (MomentumIterator it(&momentumRanges); it < DISCRETIZATION; ++it) {
			result(it.idx) = integral_screening(sc_expectation_value, it.k);
#ifndef NO_FOCK_COULOMB
			result(it.idx + DISCRETIZATION) = integral_screening(delta_n_wrapper, it.k);
#endif
#ifdef approximate_theta
			// approximate theta(omega - 0.5*|l^2 - k^2|) as theta(omega - eps_k)theta(omega - eps_l)
			if (std::abs(phonon_boundary_a(it.k) - 2. * fermi_energy) > 2.0 * omega_debye) {
				continue;
			}
#endif
			result(it.idx) -= phonon_interaction.sc_channel_integral(sc_expectation_value, it.k);
#ifndef NO_FOCK_PHONON
			result(it.idx + DISCRETIZATION) -= phonon_interaction.fock_channel_integral(delta_n_wrapper, it.k);
#endif
		}
		result(2 * DISCRETIZATION) = k_infinity_integral();

		this->Delta.fill_with(result, 0.5);
		this->Delta.clear_noise(PRECISION);
		result -= initial_values;
		++step_num;
	}

	c_float SCModel::compute_coefficient(mrock::symbolic_operators::Coefficient const& coeff, c_float first, c_float second) const
	{
		if (coeff.name == "\\epsilon_0")
		{
			return bare_dispersion(first) - fermi_energy;
		}
		else if (coeff.name == "g")
		{
#ifdef approximate_theta
			if (omega_debye >= std::abs(dispersion_to_fermi_level(first) + __fock_coulomb(first))
				&& omega_debye >= std::abs(dispersion_to_fermi_level(second) + __fock_coulomb(second)))
#else
			if (2. * omega_debye >= std::abs(phonon_boundary_a(first) - phonon_boundary_a(second)))
#endif
			{ // TODO: Think about the factor of 2, minus is handlded in the commutation program
				return phonon_coupling / rho_F;
			}
			else
			{
				return c_float{};
			}
		}
		else if (coeff.name == "V") {
			if (is_zero(coulomb_scaling)) return c_float{};
			if (coeff.momenta.front().is_zero())
				return (coulomb_scaling / PhysicalConstants::vacuum_permitivity) / (screening * screening);
			return coulomb_scaling / (2 * first * second * PhysicalConstants::vacuum_permitivity) * log_expression(first + second, first - second);
		}
		else if (coeff.name == "\\tilde{g}") {
#ifdef PHONON_SC_CHANNEL_ONLY
			return c_float{};
#else
			return 0.5 * (phonon_interaction.fock_channel(first, second) + phonon_interaction.fock_channel(second, first));
#endif
		}
		else if (coeff.name == "G") {
#ifdef PHONON_SC_CHANNEL_ONLY
			return c_float{};
#else
// TODO: Think about the factor of 2, minus is handlded in the commutation program
			return phonon_coupling / rho_F;
#endif
		}
		else
		{
			throw std::invalid_argument("Coefficient not recognized! " + coeff.name);
		}
	}

	c_float SCModel::phonon_boundary_a(const c_float k) const {
#if defined(CUT_DISPERSION_NO_COULOMB) || defined(NO_FOCK_COULOMB)
		return 2. * bare_dispersion(k);
#else
		return 2. * bare_dispersion(k) + 2. * __fock_coulomb(k);
#endif
	}

	const std::map<mrock::symbolic_operators::OperatorType, std::vector<c_complex>>& SCModel::get_expectation_values() const
	{
		if (_expecs.empty()) {
			_expecs.emplace(mrock::symbolic_operators::Number_Type, std::vector<c_complex>(DISCRETIZATION));
			_expecs.emplace(mrock::symbolic_operators::SC_Type, std::vector<c_complex>(DISCRETIZATION));
		}
		for (int k = 0; k < DISCRETIZATION; ++k) {
			_expecs.at(mrock::symbolic_operators::Number_Type)[k] = this->occupation_index(k);
			_expecs.at(mrock::symbolic_operators::SC_Type)[k] = this->sc_expectation_value_index(k);
		}

		return _expecs;
	}

	c_complex SCModel::interpolate_delta(c_float k) const {
		const int index = momentumRanges.momentum_to_floor_index(k);
		if (index >= DISCRETIZATION - 1) // Assuming Delta(k) = 0 for k -> infinity
			return (index >= DISCRETIZATION ? c_complex{} : Delta[DISCRETIZATION - 1]);
		if (index < 0) // Assuming Delta(k) = 0 for k->0
			return c_complex{};
		return mrock::utility::Numerics::interpolate_from_vector<n_interpolate>(k, momentumRanges, Delta, shifted_index(index));
	}

	c_float SCModel::interpolate_fock_correction(c_float k) const {
		const int index = momentumRanges.momentum_to_floor_index(k);
		if (index >= DISCRETIZATION - 1) // Assuming Delta(k) = 0 for k -> infinity
			return (index >= DISCRETIZATION ? c_float{} : std::real(Delta[2 * DISCRETIZATION - 1]));
		if (index < 0) // Assuming Delta(k) = const for k->0
			return c_float{};
		return std::real(mrock::utility::Numerics::interpolate_from_vector<n_interpolate>(k, momentumRanges, Delta, shifted_index(index), DISCRETIZATION));
	}

    c_float SCModel::internal_energy() const
	{
		if (!is_zero(Delta[DISCRETIZATION / 2])) {
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
			+ std::to_string(phonon_coupling) + "   omega_D="
			+ std::to_string(1e3 * omega_debye) + "meV   E_F="
			+ std::to_string(fermi_energy) + "eV   alpha=" + std::to_string(coulomb_scaling)
			+ "   lambda=" + std::to_string(screening) + "sqrt(eV)";
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
				std::string str = mrock::utility::better_to_string(number, std::chars_format::fixed);
				// Remove trailing zeroes
				str.erase(str.find_last_not_of('0') + 1, std::string::npos);
				str.erase(str.find_last_not_of('.') + 1, std::string::npos);
				return str;
			}
			};

		return "T=" + improved_string(temperature)
			+ "/coulomb_scaling=" + improved_string(coulomb_scaling)
			+ "/screening=" + improved_string(screening_ratio)
			+ "/k_F=" + improved_string(fermi_wavevector)
			+ "/g=" + improved_string(phonon_coupling)
			+ "/omega_D=" + improved_string(1e3 * omega_debye) + "/";
	}

	std::vector<c_complex> SCModel::phonon_gap() const
	{
		std::vector<c_complex> ret(DISCRETIZATION);
		for (MomentumIterator it(&momentumRanges); it < DISCRETIZATION; ++it) {
#ifdef approximate_theta
			if (std::abs(phonon_boundary_a(it.k) - 2. * fermi_energy) > 2.0 * omega_debye) {
				continue;
			}
#endif
			ret[it.idx] = -phonon_interaction.sc_channel_integral(sc_expectation_value, it.k);
		}
		return ret;
	}

	std::vector<c_complex> SCModel::coulomb_gap() const
	{
		std::vector<c_complex> ret(DISCRETIZATION);
		for (MomentumIterator it(&momentumRanges); it < DISCRETIZATION; ++it) {
			ret[it.idx] = integral_screening([this](c_float q) { return this->sc_expectation_value(q); }, it.k);
		}
		return ret;
	}

	c_complex SCModel::delta_max() const
	{
		return std::abs(*std::max_element(Delta.begin(), Delta.begin() + DISCRETIZATION,
			[](const c_complex& lhs, const c_complex& rhs) {
				return std::abs(lhs) < std::abs(rhs);
			}
		));
	}

	std::vector<c_float> SCModel::single_particle_dispersion() const
	{
		std::vector<c_float> ret(DISCRETIZATION);
		for (MomentumIterator it(&momentumRanges); it < DISCRETIZATION; ++it) {
			ret[it.idx] = bare_dispersion(it.k) - fermi_energy + __fock_coulomb(it.k);
		}
		return ret;
	}

	void SCModel::set_splines()
	{
		this->get_expectation_values();
		this->occupation = SplineContainer(_expecs[mrock::symbolic_operators::Number_Type], momentumRanges.K_MIN,
			momentumRanges.LOWER_STEP, momentumRanges.INNER_STEP, momentumRanges.UPPER_STEP,
			_OUTER_DISC, _INNER_DISC);
		this->sc_expectation_value = SplineContainer(_expecs[mrock::symbolic_operators::SC_Type], momentumRanges.K_MIN,
			momentumRanges.LOWER_STEP, momentumRanges.INNER_STEP, momentumRanges.UPPER_STEP,
			_OUTER_DISC, _INNER_DISC);
	}
}