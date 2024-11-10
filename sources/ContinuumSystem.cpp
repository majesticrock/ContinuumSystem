#include "Continuum/SCModel.hpp"
#include <Utility/Selfconsistency/IterativeSolver.hpp>
#include <Utility/OutputConvenience.hpp>
#include "Continuum/ModeHelper.hpp"
#include "Continuum/Incrementer.hpp"

#include <iomanip>
#include <omp.h>

#ifndef _NO_MPI
#include <mpi.h>
#define EXIT MPI_Finalize()
#else
#define EXIT 0
#endif

#include <filesystem>
#include <algorithm>
#include <concepts>
using namespace Continuum;

#include <SymbolicOperators/WickOperator.hpp>
#include <nlohmann/json.hpp>

const std::string BASE_FOLDER = "../../data/continuum/";

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
	for (size_t i = 0U; i < in_eV.size(); ++i) {
		ret[i] = 1e3 * std::real(in_eV[i]);
	}
	return ret;
}
template<typename number>
	requires std::floating_point<number>
std::vector<number> imag_as_meV(std::vector<std::complex<number>> const& in_eV) {
	std::vector<number> ret(in_eV.size());
	for (size_t i = 0U; i < in_eV.size(); ++i) {
		ret[i] = 1e3 * std::imag(in_eV[i]);
	}
	return ret;
}

void compute_small_U_gap() {
	Utility::InputFileReader input("params/params.config");
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
		Utility::Selfconsistency::IterativeSolver<c_complex, SCModel, ModelAttributes<c_complex>> solver(&model, &model.Delta);
		solver.compute();
		const auto buffer = model.Delta.abs().as_vector();
		gap_data[U] = { entry, *std::max_element(buffer.begin(), buffer.end()) };
	}
	Utility::saveData(gap_data, BASE_FOLDER + "test/small_U_gap.dat.gz");
}

#define RANK_RANGES(x)  const double rank_range = (std::stod(argv[3]) - init.x) / n_ranks; \
						init.x += rank * rank_range; init.recompute_dependencies();

int main(int argc, char** argv) {
	if (argc < 2) {
		std::cerr << "Invalid number of arguments: Use mpirun -n <threads> <path_to_executable> <configfile>" << std::endl;
		return -1;
	}
#ifndef _NO_MPI
	// First call MPI_Init
	MPI_Init(&argc, &argv);

	// Get my rank and the number of ranks
	int rank, n_ranks;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &n_ranks);
#else
	int rank = 0;
	int n_ranks = 1;
#endif

	Utility::InputFileReader input(argv[1]);
	Continuum::set_discretization(input.getInt("discretization_points"));

	if (false) { // compute gap in a range for small g
		compute_small_U_gap();
		return EXIT;
	}

	/*
	* Setup iterations (if asked for)
	*/
	ModelInitializer init(input);

	int n_iter = argc > 4 ? std::stoi(argv[4]) : 0;
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
			const double rank_range = 1e-3 * (std::stod(argv[3]) - init.omega_debye) / n_ranks;
			init.omega_debye += rank * rank_range;
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
	//std::cout << "Rank #" << rank << ": " << init << std::endl;
	ModeHelper modes(init);
	// We also want the last data point
	if (rank == n_ranks - 1) ++n_iter;
	for (int i = 0; i < n_iter; ++i)
	{
		/*
		* Generate setup for output
		*/
		auto delta_result = modes.getModel().Delta.real().as_vector();
		const std::string output_folder = input.getString("output_folder") + "/" + "N_k=" + std::to_string(DISCRETIZATION) + "/" + modes.getModel().to_folder();
		std::filesystem::create_directories(BASE_FOLDER + output_folder);
		auto generate_comments = [&]() {
			return nlohmann::json{
				{ "time", 				Utility::time_stamp() },
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
		Utility::saveString(jDelta.dump(4), BASE_FOLDER + output_folder + "gap.json.gz");
		std::cout << "Gap data have been saved! Delta_max = " << jDelta["Delta_max"] << std::endl;

		if (false) { // compute and save the expectation values
			auto expecs = modes.getModel().get_expectation_values();
			auto ks = modes.getModel().momentumRanges.get_k_points();

			std::vector<c_float> occupations, pairs;
			occupations.reserve(ks.size());
			pairs.reserve(ks.size());
			for (const auto& x : expecs[SymbolicOperators::Number_Type]) {
				occupations.push_back(std::real(x));
			}
			for (const auto& x : expecs[SymbolicOperators::SC_Type]) {
				pairs.push_back(std::real(x));
			}

			nlohmann::json jExpecs = {
				{"ks", ks}, {"n_k", occupations}, {"f_k", pairs}
			};
			jExpecs.merge_patch(generate_comments());
			Utility::saveString(jExpecs.dump(4), BASE_FOLDER + output_folder + "expecs.json.gz");
			std::cout << "Expectation values have been saved!" << std::endl;
		}

		if (true) {
			auto resolvents = modes.computeCollectiveModes(150);
			if (!resolvents.empty()) {
				nlohmann::json jResolvents = {
					{ "resolvents", resolvents },
					{ "continuum_boundaries", modes.continuum_boundaries() }
				};
				jResolvents.merge_patch(generate_comments());
				Utility::saveString(jResolvents.dump(4), BASE_FOLDER + output_folder + "resolvents.json.gz");
			}
		}

		if (incrementer && i < n_iter - 1) {
			incrementer->increment(init);
			modes.getModel().set_new_parameters(init);
		}
	}

	return EXIT;
}
