#include "Continuum/SCModel.hpp"
#include "../../Utility/sources/Selfconsistency/IterativeSolver.hpp"
#include "../../Utility/sources/OutputConvenience.hpp"

#include <filesystem>
using namespace Continuum;

const std::string BASE_FOLDER = "../../data/continuum/";

int main(int argc, char** argv) {
	Utility::InputFileReader input("params/test.config");
	SCModel model{ModelInitializer(input)};
	std::cout<< model.info() << std::endl;

	Utility::Selfconsistency::IterativeSolver<c_complex, SCModel, ModelAttributes<c_complex>> solver(&model, &model.Delta);
	auto result = solver.computePhases().real().as_vector();

	std::filesystem::create_directories(BASE_FOLDER + "test/");
	Utility::saveData(get_k_points(), result, BASE_FOLDER + "test/gap.dat.gz");

	return 0;
}