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
#endif
#define EXIT 0

#include <filesystem>
#include <algorithm>
using namespace Continuum;

#include <SymbolicOperators/WickOperator.hpp>
#include <nlohmann/json.hpp>

const std::string BASE_FOLDER = "../../data/continuum/";

int Continuum::DISCRETIZATION = 1000;
c_float Continuum::INV_N = 1. / Continuum::DISCRETIZATION;
int Continuum::_INNER_DISC = Continuum::DISCRETIZATION / Continuum::REL_INNER_DISCRETIZATION;
int Continuum::_OUTER_DISC = Continuum::DISCRETIZATION - Continuum::_INNER_DISC;

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

int main(int argc, char** argv) {
	if (argc < 2) {
		std::cerr << "Invalid number of arguments: Use mpirun -n <threads> <path_to_executable> <configfile>" << std::endl;
		return -1;
	}
#ifndef _NO_MPI
	// First call MPI_Init
	MPI_Init(&argc, &argv);

	// Get my rank and the number of ranks
	int rank, numberOfRanks;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numberOfRanks);
#endif

	Utility::InputFileReader input(argv[1]);
	Continuum::set_discretization(input.getInt("discretization_points"));

	if (false) { // compute gap in a range for small U
		compute_small_U_gap();
		return EXIT;
	}

	/*
	* Setup iterations (if asked for)
	*/
	ModelInitializer init(input);
	ModeHelper modes(init);
	
	const int n_iter = argc > 4 ? std::stoi(argv[4]) : 1;
	std::unique_ptr<Base_Incrementer> incrementer;
	if(argc > 4) {
		const std::string inc_type = argv[2];
		const double end_param = std::stod(argv[3]);
		if(inc_type == "T") 
			incrementer = std::make_unique<Temperature_Incrementer>((end_param - init.temperature) / n_iter);
		else if(inc_type == "g")
			incrementer = std::make_unique<PhononCoupling_Incrementer>((end_param - init.phonon_coupling) / n_iter);
		else if(inc_type == "omega_D")
			incrementer = std::make_unique<DebyeFrequency_Incrementer>((end_param - init.omega_debye) / n_iter);
		else if(inc_type == "E_F")
			incrementer = std::make_unique<FermiEnergy_Incrementer>((end_param - init.fermi_energy) / n_iter);
		else if(inc_type == "coulomb")
			incrementer = std::make_unique<CoulombScaling_Incrementer>((end_param - init.coulomb_scaling) / n_iter);
		else throw std::invalid_argument("Failed incrementer parsing. Syntax: mpirun -n <threads> <executable> <parameter_file> <incrementer_type> <end_increment> <n_increments>");
	}

	for(int i = 0; i < n_iter; ++i) 
	{
		/* 
		* Generate setup for output
		*/
		auto delta_result = modes.getModel().Delta.real().as_vector();
		const std::string output_folder = input.getString("output_folder") + "/" + modes.getModel().to_folder();
		std::filesystem::create_directories(BASE_FOLDER + output_folder);
		auto generate_comments = [&]() {
			return nlohmann::json {
				{ "time", Utility::time_stamp() },
				{ "Discretization", DISCRETIZATION },
				{ "Screening lambda", _screening },
				{ "Delta_max", std::abs(*std::max_element(delta_result.begin(), delta_result.begin() + DISCRETIZATION, 
					[](decltype(delta_result)::const_reference lhs, decltype(delta_result)::const_reference rhs){
						return std::abs(lhs) < std::abs(rhs);
					})) },
				{ "k_F", modes.getModel().fermi_wavevector },
				{ "T", modes.getModel().temperature },
				{ "g", modes.getModel().phonon_coupling },
				{ "omega_D", modes.getModel().omega_debye },
				{ "E_F", modes.getModel().fermi_energy },
				{ "Coulomb scaling", modes.getModel().coulomb_scaling }
			};
		};

		/*
		* Compute and output gap data
		*/
		nlohmann::json jDelta = {
			{"ks", modes.getModel().momentumRanges.get_k_points()},
			{"Delta_Phonon", modes.getModel().phonon_gap()},
			{"Delta_Coulomb", modes.getModel().coulomb_gap()},
			{"Delta_Fock", std::vector<double>(delta_result.begin() + DISCRETIZATION, delta_result.begin() + 2 * DISCRETIZATION)},
			{"Internal energy", modes.getModel().internal_energy()}
		};
		jDelta.merge_patch(generate_comments());
		Utility::saveString(jDelta.dump(4), BASE_FOLDER + output_folder + "gap.json.gz");
		std::cout << "Gap data have been saved!" << std::endl;

		if (false) { // compute and save the single particle energies
			std::vector<std::vector<double>> data(2, std::vector<double>(Continuum::DISCRETIZATION));
			const auto step = 2. * modes.getModel().fermi_wavevector / Continuum::DISCRETIZATION;
			for(size_t i = 0U; i < data[0].size(); ++i){
				data[0][i] = 1e-6 + i * step;
				data[1][i] = modes.getModel().energy(data[0][i]);
			}
			Utility::saveData(data, BASE_FOLDER + output_folder + "one_particle_energies.dat.dz");
		}
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

		if (false) {
			auto resolvents = modes.computeCollectiveModes(150);
			if (!resolvents.empty()) {
				nlohmann::json jResolvents = {
					{ "resolvents", resolvents },
					{ "Continuum Boundaries", modes.getModel().continuum_boundaries() }
				};
				jResolvents.merge_patch(generate_comments());
				Utility::saveString(jResolvents.dump(4), BASE_FOLDER + output_folder + "resolvents.json.gz");
			}
		}

		if (incrementer) {
			incrementer->increment(init);
			modes.getModel().set_new_parameters(init);
		}
	}

	return EXIT;
}