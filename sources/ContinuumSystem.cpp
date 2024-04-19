#include "Continuum/SCModel.hpp"
#include "../../Utility/sources/Selfconsistency/IterativeSolver.hpp"
using namespace Continuum;

int main(int argc, char** argv) {
	SCModel model({ 0, 0.1, 10 });

	Utility::Selfconsistency::IterativeSolver<c_complex, SCModel, ModelAttributes<c_complex>> solver(&model, &model.Delta);
	std::cout << solver.computePhases().real() << std::endl;

	return 0;
}