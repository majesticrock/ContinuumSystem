#include "ModelInitializer.hpp"
#include <Utility/Numerics/Roots/Bisection.hpp>
#include <Utility/ConstexprPower.hpp>

namespace Continuum {
	c_float ModelInitializer::compute_fermi_energy() const
	{
		const c_float kinetic = 0.5 * fermi_wavevector * fermi_wavevector;
		const c_float fock = (-PhysicalConstants::em_factor * SCALE) * (coulomb_scaling / r_s) *
			(1. + 0.5 * _screening * _screening * std::log(1. + 4. / (_screening * _screening)) - _screening * std::atan(2. / _screening));
		return kinetic + fock;
	}

	c_float ModelInitializer::compute_fermi_wavevector() const
	{
		return SCALE / r_s;
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