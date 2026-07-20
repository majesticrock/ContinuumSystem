import numpy as np
import matplotlib.pyplot as plt
from mrock.get_data import *
import mrock.continued_fraction as cf
import sys

data_loader = DataLoader()
    
def load_data(name):
    """ Returns the data for the Continuum test with the name 'name'.
    no_coulomb: No Coulomb interaction
    normal_screening: Normal screening
    weak_screening: Weak screening
    strong_attraction: Strong attraction
    """
    
    if name == "no_coulomb":
        return data_loader.load_panda("continuum", "test", "resolvents.json.gz",
                    **continuum_params(N_k=4000, T=0, coulomb_scaling=0, screening=0, k_F=4.25, g=0.5, omega_D=10))
    elif name == "normal_screening":
        return data_loader.load_panda("continuum", "test", "resolvents.json.gz",
                    **continuum_params(N_k=4000, T=0, coulomb_scaling=1, screening=1, k_F=4.25, g=0.8, omega_D=10))
    elif name == "weak_screening":
        return data_loader.load_panda("continuum", "test", "resolvents.json.gz",
                    **continuum_params(N_k=4000, T=0, coulomb_scaling=1, screening=1e-4, k_F=4.25, g=1, omega_D=10))
    elif name == "strong_attraction":
        return data_loader.load_panda("continuum", "test", "resolvents.json.gz",
                    **continuum_params(N_k=4000, T=0, coulomb_scaling=0, screening=0, k_F=4.25, g=1.5, omega_D=10))
    else:
        raise ValueError(f"Continuum test: Invalid name given: {name}")

def create_plot(name):
    pd_data = load_data(name)
    resolvents = cf.ContinuedFraction(pd_data, ignore_first=30, ignore_last=90)

    fig, ax = plt.subplots()
    ax.set_ylim(-0.05, 1.)
    ax.set_xlabel(r"$\omega [\mathrm{meV}]$")
    ax.set_ylabel(r"$\mathcal{A} (\omega) [\mathrm{eV}^{-1}]$")

    w_lin = np.linspace(-0.005 * pd_data["continuum_boundaries"][1], 1.1 * pd_data["continuum_boundaries"][1], 5000, dtype=complex)
    w_lin += 1e-4j

    ax.plot(1e3 * w_lin.real, resolvents.spectral_density(w_lin, "phase_SC",     with_terminator=True), label="Phase")
    ax.plot(1e3 * w_lin.real, resolvents.spectral_density(w_lin, "amplitude_SC", with_terminator=True), label="Higgs")

    resolvents.mark_continuum(ax, 1e3)

    ax.set_xlim(1e3 * np.min(w_lin.real), 1e3 * np.max(w_lin.real))
    ax.set_title(f"Continuum: {name}")
    ax.legend()
    fig.tight_layout()
    plt.show()
    
if len(sys.argv) > 1:
    create_plot(sys.argv[1])
else:
    print("Please provide the name of the test you would like to plot.")