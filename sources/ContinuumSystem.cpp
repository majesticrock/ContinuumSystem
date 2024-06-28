#include "Continuum/SCModel.hpp"
#include <Utility/Selfconsistency/IterativeSolver.hpp>
#include <Utility/OutputConvenience.hpp>
#include "Continuum/ModeHelper.hpp"

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

const std::string BASE_FOLDER = "../../data/continuum/";

int Continuum::DISCRETIZATION = 1000;
c_float Continuum::INV_N = 1. / Continuum::DISCRETIZATION;
int Continuum::_INNER_DISC = Continuum::DISCRETIZATION / Continuum::REL_INNER_DISCRETIZATION;
int Continuum::_OUTER_DISC = Continuum::DISCRETIZATION - Continuum::_INNER_DISC;


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
		return EXIT;
	}

	ModeHelper modes(input);
	auto delta_result = modes.getModel().Delta.real().as_vector();

	std::vector<std::string> comments;
	comments.push_back("Discretization: " + std::to_string(DISCRETIZATION));
	comments.push_back("k_F: " + std::to_string(modes.getModel().fermi_wavevector));

	Utility::saveData(std::vector<std::vector<double>>{
		modes.getModel().momentumRanges.get_k_points(),
			modes.getModel().phonon_gap(),
			modes.getModel().coulomb_gap(),
			std::vector<double>(delta_result.begin() + DISCRETIZATION, delta_result.begin() + 2 * DISCRETIZATION)
	}, BASE_FOLDER + output_folder + "/gap.dat.gz", comments);
	std::cout << "Gap data have been saved!" << std::endl;

	std::cout << "Delta_max = " << std::scientific << std::setprecision(14)
		<< *std::max_element(delta_result.begin(), delta_result.begin() + DISCRETIZATION) << std::endl;
	std::cout << "Delta_min = " << std::scientific << std::setprecision(14)
		<< *std::min_element(delta_result.begin(), delta_result.begin() + DISCRETIZATION) << std::endl;
	std::cout << "Internal energy = " << modes.getModel().internal_energy() << std::endl;

	Utility::saveData(modes.getModel().continuum_boundaries(), BASE_FOLDER + output_folder + "/continuum.dat.gz");

	if constexpr (true) { // compute and save the single particle energies
		std::vector<std::vector<double>> data(2, std::vector<double>(Continuum::DISCRETIZATION));
		const auto step = 2. * modes.getModel().fermi_wavevector / Continuum::DISCRETIZATION;
		for(size_t i = 0U; i < data[0].size(); ++i){
			data[0][i] = 1e-6 + i * step;
			data[1][i] = modes.getModel().energy(data[0][i]);
		}
		Utility::saveData(data, BASE_FOLDER + output_folder + "/one_particle_energies.dat.dz");
	}

	if constexpr (true) { // compute and save the expectation values
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

		Utility::saveData(std::vector<std::vector<c_float>>{ks, occupations, pairs}, BASE_FOLDER + output_folder + "/expecs.dat.gz");
		std::cout << "Expectation values have been saved!" << std::endl;
	}

	return EXIT;
	auto mode_result = modes.computeCollectiveModes(150);
	if (!mode_result.empty()) {
#ifdef _complex
		std::vector<std::string> names { 
					"phase_SC_a", "phase_SC_a+b", "phase_SC_a+ib",
					"higgs_SC_a", "higgs_SC_a+b", "higgs_SC_a+ib"
					};
#else
		std::vector<std::string> names{ "phase_SC", "higgs_SC" };
#endif
		for (size_t i = 0U; i < mode_result.size(); ++i)
		{
			mode_result[i].writeDataToFile(BASE_FOLDER + output_folder + "/resolvent_" + names[i], comments);
		}
	}

	return EXIT;
}