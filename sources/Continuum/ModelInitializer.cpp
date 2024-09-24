#include "ModelInitializer.hpp"
#include <Utility/Numerics/Roots/Bisection.hpp>
#include <Utility/ConstexprPower.hpp>

namespace Continuum {
	c_float ModelInitializer::compute_screening() const
	{
		return screening_ratio * PhysicalConstants::screening_prefactor * sqrt(fermi_wavevector);
	}
	c_float ModelInitializer::compute_fermi_energy() const
	{
		const c_float kinetic = bare_dispersion(fermi_wavevector);
		if(is_zero(coulomb_scaling)) return kinetic;
		const c_float lambda = screening / fermi_wavevector;
		const c_float fock = - PhysicalConstants::em_factor * coulomb_scaling * fermi_wavevector * 
			(1. + 0.5 * lambda * lambda * std::log(1. + 4. / (lambda * lambda)) - lambda * std::atan(2. / lambda));
		return kinetic + fock;
	}

	c_float ModelInitializer::compute_rho_F() const
	{
		return fermi_wavevector / (2. * PI * PI);
		/* const c_float kf2 = fermi_wavevector * fermi_wavevector;
		const c_float ks2 = screening * screening;
		return kf2 / (2. * PI * PI * (fermi_wavevector - 0.5 * PhysicalConstants::em_factor * coulomb_scaling * 
			( 2. * kf2 - ks2 ) * ( std::log(1. + 4. * kf2 / ks2) / kf2 - 4. / (4. * kf2 + ks2) )
		));*/
	}

	std::ostream& operator<<(std::ostream& os, ModelInitializer const& init) 
	{
		os << "T=" << init.temperature << " K   "
		<< "g=" << init.phonon_coupling << "   "
		<< "omega_D=" << init.omega_debye << " eV   "
		<< "E_F=" << init.fermi_energy << " eV   "
		<< "alpha=" << init.coulomb_scaling << "   "
		<< "k_F=" << init.fermi_wavevector << " sqrt(eV)   "
		<< "rho_F=" << init.rho_F << " sqrt(eV)   ";
	return os;
}
}