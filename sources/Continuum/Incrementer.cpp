#include "Incrementer.hpp"

namespace Continuum {
    void Temperature_Incrementer::increment(ModelInitializer& init) const 
        { init.temperature += _Delta; }
    double Temperature_Incrementer::current(ModelInitializer const& init) const
        { return init.temperature; }

    void PhononCoupling_Incrementer::increment(ModelInitializer& init) const 
        { init.phonon_coupling += _Delta; }
    double PhononCoupling_Incrementer::current(ModelInitializer const& init) const
        { return init.phonon_coupling; }

    void DebyeFrequency_Incrementer::increment(ModelInitializer& init) const
        { init.omega_debye += _Delta; }
    double DebyeFrequency_Incrementer::current(ModelInitializer const& init) const
        { return init.omega_debye; }

    void FermiEnergy_Incrementer::increment(ModelInitializer& init) const
        { init.fermi_energy += _Delta; init.compute_fermi_wavevector(); init.compute_v_over_n(); }
    double FermiEnergy_Incrementer::current(ModelInitializer const& init) const
        { return init.fermi_energy; }

    void CoulombScaling_Incrementer::increment(ModelInitializer& init) const
        { init.coulomb_scaling += _Delta; init.compute_fermi_wavevector(); init.compute_v_over_n(); }
    double CoulombScaling_Incrementer::current(ModelInitializer const& init) const
        { return init.coulomb_scaling; }
}