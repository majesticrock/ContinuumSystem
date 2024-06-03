#include "Continuum/SCModel.hpp"
#include "../../Utility/sources/Selfconsistency/IterativeSolver.hpp"
#include "../../Utility/sources/OutputConvenience.hpp"
#include "Continuum/ModeHelper.hpp"

#include <iomanip>
#include <omp.h>

#define _NO_MPI
#ifndef _NO_MPI
#include <mpi.h>
#define EXIT MPI_Finalize()
#endif
#define EXIT 0

#include <filesystem>
#include <algorithm>
using namespace Continuum;

#include "../../FermionCommute/sources/WickOperator.hpp"

const std::string BASE_FOLDER = "../../data/continuum/";

int Continuum::DISCRETIZATION = 500;
c_float Continuum::INV_N = 1. / Continuum::DISCRETIZATION;

int main(int argc, char** argv) {
#ifndef _NO_MPI
	// First call MPI_Init
	MPI_Init(&argc, &argv);

	// Get my rank and the number of ranks
	int rank, numberOfRanks;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numberOfRanks);
#endif

	Utility::InputFileReader input("params/params.config");
	const std::string output_folder = input.getString("output_folder");
	Continuum::set_discretization(input.getInt("discretization_points"));
	std::filesystem::create_directories(BASE_FOLDER + output_folder + "/");

	if constexpr (false) { // compute gap in a range for small U
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
			//std::cout << model.info() << "\t" << *std::max_element(buffer.begin(), buffer.end()) << std::endl;
		}
		Utility::saveData(gap_data, BASE_FOLDER + "test/small_U_gap.dat.gz");
		return 0;
	}

	ModeHelper modes(input);
	auto delta_result = modes.getModel().Delta.real().as_vector();
	auto delta_imag = modes.getModel().Delta.imag().as_vector();
	//Utility::saveData(modes.getModel().get_k_points(), delta_result, BASE_FOLDER + output_folder + "/gap.dat.gz");
	Utility::saveData(std::vector<std::vector<double>>{
		modes.getModel().get_k_points(),
			std::vector<double>(delta_result.begin(), delta_result.begin() + DISCRETIZATION),
			std::vector<double>(delta_imag.begin(), delta_imag.begin() + DISCRETIZATION)
			//std::vector<double>(delta_result.begin() + DISCRETIZATION, delta_result.begin() + 2 * DISCRETIZATION),
	}, BASE_FOLDER + output_folder + "/gap.dat.gz");
	std::cout << "Gap data have been saved!" << std::endl;

	std::cout << "Delta_max = " << std::scientific << std::setprecision(14)
		<< *std::max_element(delta_result.begin(), delta_result.begin() + DISCRETIZATION) << std::endl;
	std::cout << "Delta_min = " << std::scientific << std::setprecision(14)
		<< *std::min_element(delta_result.begin(), delta_result.begin() + DISCRETIZATION) << std::endl;
	std::cout << "Internal energy = " << modes.getModel().internal_energy() << std::endl;

	Utility::saveData(modes.getModel().continuum_boundaries(), BASE_FOLDER + output_folder + "/continuum.dat.gz");

	if constexpr (true) { // compute and save the expectation values
		auto expecs = modes.getModel().get_expectation_values();
		auto ks = modes.getModel().get_k_points();

		std::vector<c_float> occupations, pairs;
		occupations.reserve(ks.size());
		pairs.reserve(ks.size());
		for (const auto& x : expecs[SymbolicOperators::Number_Type]) {
			occupations.push_back(std::real(x));
		}
		for (const auto& x : expecs[SymbolicOperators::SC_Type]) {
			pairs.push_back(std::real(x));
		}

		Utility::saveData(std::vector<std::vector<c_float>>{ks, occupations, pairs}, BASE_FOLDER + output_folder + "/expecs.dat.gz");
		std::cout << "Expectation values have been saved!" << std::endl;
	}
	return EXIT;

	auto mode_result = modes.computeCollectiveModes(150);
	if (!mode_result.empty()) {
		std::vector<std::string> comments;
		comments.push_back("Discretization: " + std::to_string(DISCRETIZATION));

		//std::vector<std::string> names{ "higgs_SC_a", "higgs_SC_a+b", "higgs_SC_a+ib",
		//			"phase_SC_a", "phase_SC_a+b", "phase_SC_a+ib" };
		std::vector<std::string> names{ "phase_SC", "higgs_SC" };

		for (size_t i = 0U; i < mode_result.size(); ++i)
		{
			mode_result[i].writeDataToFile(BASE_FOLDER + output_folder + "/resolvent_" + names[i], comments);
		}
	}

	return EXIT;
}