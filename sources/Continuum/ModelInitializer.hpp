#pragma once
#include "GlobalDefinitions.hpp"
#include <Utility/InputFileReader.hpp>

namespace Continuum {
	struct ModelInitializer {
		// Members set by the parameter file
		c_float temperature;
		c_float phonon_coupling;
		c_float omega_debye;
		c_float fermi_wavevector;
		c_float coulomb_scaling;
		c_float screening_ratio;

		// Members computed depending on set members
		c_float screening;
		c_float fermi_energy;
		c_float rho_F;

		ModelInitializer(Utility::InputFileReader& input)
			: temperature{ input.getDouble("T") },
			phonon_coupling{ input.getDouble("phonon_coupling") },
			omega_debye{ 1e-3 * input.getDouble("omega_debye") }, // given in meV in the parameter file
			fermi_wavevector{ input.getDouble("k_F") },
			coulomb_scaling{ input.getDouble("coulomb_scaling") },
			screening_ratio{ is_zero(coulomb_scaling) ? c_float{} : input.getDouble("screening") },
			screening{ compute_screening() },
			fermi_energy{ compute_fermi_energy() },
			rho_F{ compute_rho_F() }
		{ };

		c_float compute_screening() const;
		c_float compute_fermi_energy() const;
		c_float compute_rho_F() const;
		inline void recompute_dependencies() {
			this->screening = compute_screening();
			this->fermi_energy = compute_fermi_energy();
			this->rho_F = compute_rho_F();
		}
	};

	std::ostream& operator<<(std::ostream& os, ModelInitializer const& init);
}