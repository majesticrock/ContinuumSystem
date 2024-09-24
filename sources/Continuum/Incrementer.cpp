#include "Incrementer.hpp"

namespace Continuum {
	void Temperature_Incrementer::increment(ModelInitializer& init, int n/*=1*/) const
	{
		init.temperature = n * _Delta + init.temperature;
	}
	double Temperature_Incrementer::current(ModelInitializer const& init) const
	{
		return init.temperature;
	}

	void PhononCoupling_Incrementer::increment(ModelInitializer& init, int n/*=1*/) const
	{
		init.phonon_coupling = n * _Delta + init.phonon_coupling;
	}
	double PhononCoupling_Incrementer::current(ModelInitializer const& init) const
	{
		return init.phonon_coupling;
	}

	void DebyeFrequency_Incrementer::increment(ModelInitializer& init, int n/*=1*/) const
	{
		init.omega_debye = n * _Delta + init.omega_debye;
	}
	double DebyeFrequency_Incrementer::current(ModelInitializer const& init) const
	{
		return init.omega_debye;
	}

	void FermiWavevector_Incrementer::increment(ModelInitializer& init, int n/*=1*/) const
	{
		init.fermi_wavevector = n * _Delta + init.fermi_wavevector;
		init.fermi_energy = init.compute_fermi_energy();
	}
	double FermiWavevector_Incrementer::current(ModelInitializer const& init) const
	{
		return init.fermi_wavevector;
	}

	void CoulombScaling_Incrementer::increment(ModelInitializer& init, int n/*=1*/) const
	{
		init.coulomb_scaling = n * _Delta + init.coulomb_scaling;
		init.fermi_energy = init.compute_fermi_energy();
	}
	double CoulombScaling_Incrementer::current(ModelInitializer const& init) const
	{
		return init.coulomb_scaling;
	}
}