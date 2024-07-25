#pragma once
#include "GlobalDefinitions.hpp"
#include <Utility/InputFileReader.hpp>

namespace Continuum {
	struct ModelInitializer {
		// Members set by the parameter file
		c_float temperature;
		c_float phonon_coupling;
		c_float omega_debye;
		c_float fermi_energy;
		c_float coulomb_scaling;

		// Members computed depending on set members
		c_float fermi_wavevector;

		ModelInitializer(Utility::InputFileReader& input)
			: temperature{ PhysicalConstants::k_B * input.getDouble("T") }, 
			phonon_coupling{ input.getDouble("phonon_coupling") },
			omega_debye{ 1e-3 * input.getDouble("omega_debye") }, // given in meV in the parameter file
			fermi_energy{ input.getDouble("fermi_energy") },
			coulomb_scaling{ input.getDouble("coulomb_scaling") },
			fermi_wavevector{ compute_fermi_wavevector() }
		{ };

		c_float compute_fermi_wavevector() const;
		inline void recompute_dependencies() {
			this->fermi_wavevector = compute_fermi_wavevector();
		}
	};

	std::ostream& operator<<(std::ostream& os, ModelInitializer const& init);
}