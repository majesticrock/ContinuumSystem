#include "Continuum/SCModel.hpp"
#include "../../Utility/sources/Selfconsistency/IterativeSolver.hpp"
#include "../../Utility/sources/OutputConvenience.hpp"
#include "Continuum/ModeHelper.hpp"

#include <filesystem>
using namespace Continuum;

const std::string BASE_FOLDER = "../../data/continuum/";

int main(int argc, char** argv) {
	Utility::InputFileReader input("params/test.config");

	ModeHelper modes(input);

	std::cout << modes.getModel().info() << std::endl;
	auto delta_result = modes.getModel().Delta.real().as_vector();
	std::filesystem::create_directories(BASE_FOLDER + "test/");
	Utility::saveData(get_k_points(), delta_result, BASE_FOLDER + "test/gap.dat.gz");

	auto mode_result = modes.computeCollectiveModes(150);
	if (!mode_result.empty()) {
		std::vector<std::string> comments;
		comments.push_back("Discretization: " + std::to_string(DISCRETIZATION));

		std::vector<std::string> names{ "higgs_SC_a", "higgs_SC_a+b", "higgs_SC_a+ib",
					"phase_SC_a", "phase_SC_a+b", "phase_SC_a+ib" };

		for (size_t i = 0U; i < mode_result.size(); ++i)
		{
			mode_result[i].writeDataToFile(BASE_FOLDER + "test/resolvent_" + names[i], comments);
		}
	}

	return 0;
}