#include "Continuum/SCModel.hpp"
#include "../../Utility/sources/Selfconsistency/IterativeSolver.hpp"
#include "../../Utility/sources/OutputConvenience.hpp"
#include "Continuum/ModeHelper.hpp"

#include <omp.h>

#include <filesystem>
#include <algorithm>
using namespace Continuum;

#include "../../FermionCommute/sources/WickOperator.hpp"

const std::string BASE_FOLDER = "../../data/continuum/";

int Continuum::DISCRETIZATION = 500;
c_float Continuum::INV_N = 1. / Continuum::DISCRETIZATION;

int main(int argc, char** argv) {
	Utility::InputFileReader input("params/test.config");
	Continuum::set_discretization(input.getInt("discretization_points"));
	std::filesystem::create_directories(BASE_FOLDER + "test/");

	if (false) { // compute gap in a range for small U
		constexpr int N_points = 200;
		constexpr double step = 0.5;
		constexpr double begin = 150;
		std::vector<std::vector<c_float>> gap_data(N_points);
		
#pragma omp parallel for
		for (int U = 0; U < N_points; ++U) {
			ModelInitializer init(input);
			double entry = begin + step * U;
			init.U = entry;
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

	auto delta_result = modes.getModel().Delta.abs().as_vector();
	Utility::saveData(modes.getModel().get_k_points(), delta_result, BASE_FOLDER + "test/gap.dat.gz");
	std::cout << "Gap data have been saved!" << std::endl;
	std::cout << "Delta_max = " << *std::max_element(delta_result.begin(), delta_result.end()) << std::endl;

	Utility::saveData(modes.getModel().continuum_boundaries(), BASE_FOLDER + "test/continuum.dat.gz");

	if (true) {
		auto expecs = modes.getModel().get_expectation_values();
		auto ks = modes.getModel().get_k_points();

		std::vector<double> occupations, pairs;
		occupations.reserve(ks.size());
		pairs.reserve(ks.size());
		for(const auto& x : expecs[SymbolicOperators::Number_Type]){
			occupations.push_back(x.real());
		}
		for(const auto& x : expecs[SymbolicOperators::SC_Type]){
			pairs.push_back(x.real());
		}

		Utility::saveData(std::vector<std::vector<double>>{ks, occupations, pairs}, BASE_FOLDER + "test/expecs.dat.gz");
	}

	auto mode_result = modes.computeCollectiveModes(150);
	if (!mode_result.empty()) {
		std::vector<std::string> comments;
		comments.push_back("Discretization: " + std::to_string(DISCRETIZATION));

		//std::vector<std::string> names{ "higgs_SC_a", "higgs_SC_a+b", "higgs_SC_a+ib",
		//			"phase_SC_a", "phase_SC_a+b", "phase_SC_a+ib" };
		std::vector<std::string> names{ "higgs_SC", "phase_SC" };

		for (size_t i = 0U; i < mode_result.size(); ++i)
		{
			mode_result[i].writeDataToFile(BASE_FOLDER + "test/resolvent_" + names[i], comments);
		}
	}

	return 0;
}