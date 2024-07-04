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
		c_float V_OVER_N;

		ModelInitializer(Utility::InputFileReader& input)
			: temperature{ PhysicalConstants::k_B * input.getDouble("T") }, phonon_coupling{ input.getDouble("phonon_coupling") },
			omega_debye{ input.getDouble("omega_debye") }, fermi_energy{ input.getDouble("fermi_energy") },
			coulomb_scaling{ input.getDouble("coulomb_scaling") },
			fermi_wavevector{ compute_fermi_wavevector() },
			V_OVER_N{ compute_v_over_n() }
		{ };

		c_float compute_fermi_wavevector() const;
		c_float compute_v_over_n() const;
	};
}