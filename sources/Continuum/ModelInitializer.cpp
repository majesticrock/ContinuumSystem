#include "ModelInitializer.hpp"
#include <Utility/Numerics/Roots/Bisection.hpp>
#include <Utility/ConstexprPower.hpp>

namespace Continuum {
	c_float ModelInitializer::compute_fermi_wavevector() const
	{
		c_float kF = coulomb_scaling * PhysicalConstants::em_factor 
			+ sqrt((coulomb_scaling * PhysicalConstants::em_factor * coulomb_scaling * PhysicalConstants::em_factor) + 2. * PhysicalConstants::effective_mass * fermi_energy);

		auto energy_at_fermi_level = [this](c_float kF) {
			const c_float k_sum{ 2 * kF };
			const c_float ln_factor{ (_screening * _screening) / (2.0 * kF * kF) };
			const c_float fock = -coulomb_scaling * PhysicalConstants::em_factor * kF * 
			(
				1.0 + ln_factor * log_expression(k_sum, c_float{}) - (_screening / kF) * std::atan(k_sum / _screening)
			);
			return bare_dispersion(kF) + fock - fermi_energy;
		};
		return Utility::Numerics::Roots::bisection(energy_at_fermi_level, 0.95 * kF, 1.05 * kF, PRECISION, 250);
	}

	std::ostream& operator<<(std::ostream& os, ModelInitializer const& init) 
	{
		os << "T=" << init.temperature << " K   "
		<< "g=" << init.phonon_coupling << " eV   "
		<< "omega_D=" << init.omega_debye << " eV   "
		<< "E_F=" << init.fermi_energy << " eV   "
		<< "alpha=" << init.coulomb_scaling << "   "
		<< "k_F=" << init.fermi_wavevector << " sqrt(eV)   ";
	return os;
}
}