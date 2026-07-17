#include <mrock/utility/Selfconsistency/IterativeSolver.hpp>
#include <mrock/utility/OutputConvenience.hpp>
#include <mrock/utility/info_to_json.hpp>
#include <mrock/symbolic_operators/WickOperator.hpp>
#include <nlohmann/json.hpp>
#include <mrock/info.h>

#include <iomanip>
#include <omp.h>
#include <filesystem>
#include <algorithm>
#include <concepts>

#include "Continuum/SCModel.hpp"
#include "Continuum/ModeHelper.hpp"
#include "Continuum/Incrementer.hpp"
// File is generated on build by cmake
#include <continuum/info.h>

#ifndef OUTPUT_DATA_DIR
#define OUTPUT_DATA_DIR "../../data/continuum/"
#endif

using namespace Continuum;

int Continuum::DISCRETIZATION = 1000;
c_float Continuum::INV_N = 1. / Continuum::DISCRETIZATION;
int Continuum::_INNER_DISC = Continuum::DISCRETIZATION / Continuum::REL_INNER_DISCRETIZATION;
int Continuum::_OUTER_DISC = Continuum::DISCRETIZATION - Continuum::_INNER_DISC;

template<typename number>
	requires std::floating_point<number>
constexpr number as_meV(number in_eV) {
	in_eV *= 1e3;
	return in_eV;
}
template<typename number>
	requires std::floating_point<number>
std::vector<number>&& as_meV(std::vector<number>&& in_eV) {
	std::ranges::for_each(in_eV, [](number& num) { num *= 1e3; });
	return std::move(in_eV);
}
template<typename number>
	requires std::floating_point<number>
std::vector<number> as_meV(std::vector<std::complex<number>> const& in_eV) {
	std::vector<number> ret(in_eV.size());
	for (std::size_t i = 0U; i < in_eV.size(); ++i) {
		ret[i] = 1e3 * std::real(in_eV[i]);
	}
	return ret;
}
template<typename number>
	requires std::floating_point<number>
std::vector<number> imag_as_meV(std::vector<std::complex<number>> const& in_eV) {
	std::vector<number> ret(in_eV.size());
	for (std::size_t i = 0U; i < in_eV.size(); ++i) {
		ret[i] = 1e3 * std::imag(in_eV[i]);
	}
	return ret;
}

void compute_small_U_gap() {
	mrock::utility::InputFileReader input("params/params.config");
	constexpr int N_points = 200;
	constexpr c_float step = 0.03;
	constexpr c_float begin = 1;
	std::vector<std::vector<c_float>> gap_data(N_points);

#pragma omp parallel for
	for (int U = 0; U < N_points; ++U) {
		ModelInitializer init(input);
		c_float entry = begin + step * U;
		init.phonon_coupling = entry;
		SCModel model(init);
		mrock::utility::Selfconsistency::IterativeSolver<c_complex, SCModel, ModelAttributes<c_complex>> solver(&model, &model.Delta);
		solver.compute();
		const auto buffer = model.Delta.abs().as_vector();
		gap_data[U] = { entry, *std::max_element(buffer.begin(), buffer.end()) };
	}
	mrock::utility::save_data(gap_data, std::string(OUTPUT_DATA_DIR) + "test/small_U_gap.dat.gz");
}

#define RANK_RANGES(x)  const double rank_range = (std::stod(argv[3]) - init.x); \
						init.recompute_dependencies();

int main(int argc, char** argv) {
	if (argc < 2) {
		std::cerr << "Invalid number of arguments: Use ./path/to/executable <configfile>" << std::endl;
		return -1;
	}

	mrock::utility::InputFileReader input(argv[1]);
	Continuum::set_discretization(input.getInt("discretization_points"));

	if (false) { // compute gap in a range for small g
		compute_small_U_gap();
		return 0;
	}

	/*
	* Setup iterations (if asked for)
	*/
	ModelInitializer init(input);

	// +1 to also include the last data point
	int n_iter = (argc > 4 ? std::stoi(argv[4]) : 0) + 1;
	std::unique_ptr<Base_Incrementer> incrementer;
	if (argc > 4) {
		const std::string inc_type = argv[2];
		if (inc_type == "T" || inc_type == "temperature")
		{
			RANK_RANGES(temperature);
			incrementer = std::make_unique<Temperature_Incrementer>(rank_range / n_iter);
		}
		else if (inc_type == "g" || inc_type == "phonon_coupling")
		{
			RANK_RANGES(phonon_coupling);
			incrementer = std::make_unique<PhononCoupling_Incrementer>(rank_range / n_iter);
		}
		else if (inc_type == "omega_D" || inc_type == "omega_debye")
		{
			const double rank_range = 1e-3 * (std::stod(argv[3]) - init.omega_debye);
			incrementer = std::make_unique<DebyeFrequency_Incrementer>(rank_range / n_iter);
		}
		else if (inc_type == "k_F" || inc_type == "fermi_wavevector")
		{
			RANK_RANGES(fermi_wavevector);
			incrementer = std::make_unique<FermiWavevector_Incrementer>(rank_range / n_iter);
		}
		else if (inc_type == "coulomb" || inc_type == "coulomb_scaling")
		{
			RANK_RANGES(coulomb_scaling);
			incrementer = std::make_unique<CoulombScaling_Incrementer>(rank_range / n_iter);
		}
		else throw std::invalid_argument("Failed incrementer parsing. Syntax: mpirun -n <threads> <executable> <parameter_file> <incrementer_type> <end_increment> <n_increments>");
	}

	ModeHelper modes(init);
	for (int i = 0; i < n_iter; ++i)
	{
		/*
		* Generate setup for output
		*/
		auto delta_result = modes.getModel().Delta.real().as_vector();
		const std::string output_folder = input.getString("output_folder") + "/" + "N_k=" + std::to_string(DISCRETIZATION) + "/" + modes.getModel().to_folder();
		std::filesystem::create_directories(std::string(OUTPUT_DATA_DIR) + output_folder);
		auto generate_comments = [&]() {
			return nlohmann::json{
				{ "time", 				mrock::utility::time_stamp() },
				{ "discretization", 	DISCRETIZATION },
				{ "inner_discretization", _INNER_DISC },
				{ "lambda_screening", 	modes.getModel().screening_ratio },
				{ "Delta_max", 			as_meV(modes.getModel().delta_max()) },
				{ "k_F", 				modes.getModel().fermi_wavevector },
				{ "T", 					modes.getModel().temperature },
				{ "g", 					modes.getModel().phonon_coupling },
				{ "omega_D", 			as_meV(modes.getModel().omega_debye) },
				{ "E_F", 				modes.getModel().fermi_energy },
				{ "coulomb_scaling",	modes.getModel().coulomb_scaling },
				{ "k_infinity_factor", 	std::real(2. * PhysicalConstants::em_factor * modes.getModel().coulomb_scaling * delta_result[2 * DISCRETIZATION]) },
				{ "k_zero_factor", 		std::real(modes.getModel().k_zero_integral()) },
				{ "internal_energy", 	modes.getModel().internal_energy() }
			};
			};

		// Generate metadata
		nlohmann::json info_json = mrock::utility::generate_json<ContinuumSystem::info>("continuum_");
		info_json.update(mrock::utility::generate_json<mrock::info>("mrock_"));
		mrock::utility::save_string(info_json.dump(4), std::string(OUTPUT_DATA_DIR) + output_folder + "metadata.json.gz");

		/*
		* Compute and output gap data
		*/
		nlohmann::json jDelta = generate_comments();
		jDelta.update(nlohmann::json{
			{ "data", {
#ifdef _complex
					{ "imag_Delta_Phonon", 		imag_as_meV(modes.getModel().phonon_gap()) },
					{ "imag_Delta_Coulomb", 	imag_as_meV(modes.getModel().coulomb_gap()) },
#endif
					{ "ks", 			modes.getModel().momentumRanges.get_k_points() },
					{ "Delta_Phonon", 	as_meV(modes.getModel().phonon_gap()) },
					{ "Delta_Coulomb", 	as_meV(modes.getModel().coulomb_gap()) },
					{ "Delta_Fock", 	as_meV(std::vector<double>(delta_result.begin() + DISCRETIZATION, delta_result.begin() + 2 * DISCRETIZATION)) },
					{ "xis", 			modes.getModel().single_particle_dispersion() }
				}
			}
			});
		mrock::utility::save_string(jDelta.dump(4), std::string(OUTPUT_DATA_DIR) + output_folder + "gap.json.gz");
		std::cout << "Gap data have been saved! Delta_max = " << jDelta["Delta_max"] << std::endl;

#ifndef CONTINUUM_FULL_DIAG
		{
			auto resolvents = modes.compute_collective_modes(150);
			if (!resolvents.empty()) {
				nlohmann::json jResolvents = {
					{ "resolvents", resolvents },
					{ "continuum_boundaries", modes.continuum_boundaries() }
				};
				jResolvents.merge_patch(generate_comments());
				mrock::utility::save_string(jResolvents.dump(4), std::string(OUTPUT_DATA_DIR) + output_folder + "resolvents.json.gz");
			}
		}
#else // CONTINUUM_FULL_DIAG
		{
			auto full_diag_data = modes.full_diagonalization();
			nlohmann::json jFullDiag = {
				{ "phase", full_diag_data.first },
				{ "amplitude", full_diag_data.second },
				{ "continuum_boundaries", modes.continuum_boundaries() }
			};
			jFullDiag.merge_patch(generate_comments());
			mrock::utility::save_string(jFullDiag.dump(4), std::string(OUTPUT_DATA_DIR) + output_folder + "full_diagonalization.json.gz");
		}
#endif // CONTINUUM_FULL_DIAG

		if (incrementer && i < n_iter - 1) {
			incrementer->increment(init);
			modes.getModel().set_new_parameters(init);
		}
	}

	return 0;
}
