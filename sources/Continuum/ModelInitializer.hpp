#pragma once
#include "GlobalDefinitions.hpp"
#include <Utility/InputFileReader.hpp>

namespace Continuum {
	struct ModelInitializer {
		// Members set by the parameter file
		c_float temperature{};
		c_float phonon_coupling{};
		c_float omega_debye{};
		c_float coulomb_scaling{};
		c_float r_s{};

		// Members computed depending on set members
		c_float fermi_wavevector;
		c_float fermi_energy;

		ModelInitializer(Utility::InputFileReader& input)
			: temperature{ PhysicalConstants::k_B * input.getDouble("T") }, 
			phonon_coupling{ input.getDouble("phonon_coupling") },
			omega_debye{ 1e-3 * input.getDouble("omega_debye") }, // given in meV in the parameter file
			coulomb_scaling{ input.getDouble("coulomb_scaling") },
			r_s{ input.getDouble("r_s") },
			fermi_wavevector{ compute_fermi_wavevector() },
			fermi_energy{ compute_fermi_energy() }
		{ };

		c_float compute_fermi_energy() const;
		c_float compute_fermi_wavevector() const;
		inline void recompute_dependencies() {
			this->fermi_wavevector = compute_fermi_wavevector();
			this->fermi_energy = compute_fermi_energy();
		}
	};

	std::ostream& operator<<(std::ostream& os, ModelInitializer const& init);
}