import numpy as np
from reference_file_tddft import *

def test_groundstate(converged_results: dict):
    """Test results in main output file INFO.OUT
    containing the essential information on the material system, parameters of the calculation and results
    """
    total_energy = -5314.76550192
    assert np.isclose(converged_results['Total energy'], total_energy),\
        f"Incorrect value for total energy. Expect {total_energy}"

    fermi_energy = 0.33114209
    assert np.isclose(converged_results['Fermi energy'], fermi_energy),\
        f"Incorrect value for fermi energy. Expect {fermi_energy}"

    kinetic_energy = 5560.16937117
    assert np.isclose(converged_results['Kinetic energy'], kinetic_energy),\
        f"Incorrect value for kinetic energy. Expect {kinetic_energy}"

    coulomb_energy = -10728.00832093
    assert np.isclose(converged_results['Coulomb energy'], coulomb_energy),\
        f"Incorrect value for coulomb energy. Expect {coulomb_energy}"

    exchange_energy = -144.07719204
    assert np.isclose(converged_results['Exchange energy'], exchange_energy),\
        f"Incorrect value for exchange energy. Expect {exchange_energy}"

    correlation_energy = -2.84936011
    assert np.isclose(converged_results['Correlation energy'], correlation_energy),\
        f"Incorrect value for correlation energy. Expect {correlation_energy}"

    fermi_energy_dos = 31.30121651
    assert np.isclose(converged_results['DOS at Fermi energy (states/Ha/cell)'], fermi_energy_dos),\
        f"Incorrect value for DOS at Fermi energy (states/Ha/cell). Expect {fermi_energy_dos}"

    charge_core_electrons = 28.0
    assert np.isclose(converged_results['core'], charge_core_electrons),\
        f"Incorrect value for charge of core electrons. Expect {charge_core_electrons}"

    core_leakage_charge = 6.2e-07
    assert np.isclose(converged_results['core leakage'],core_leakage_charge),\
        f"Incorrect value for core leakage charge. Expect {core_leakage_charge}"

    charge_valence_electrons = 19.0
    assert np.isclose(converged_results['valence'], charge_valence_electrons),\
        f"Incorrect value for charge of valence electrons. Expect {charge_valence_electrons}"

    interstitial_region_charge = 2.72325304
    assert np.isclose(converged_results['interstitial'], interstitial_region_charge),\
        f"Incorrect value for charge in interstitial region. Expect {interstitial_region_charge}"

    muffin_tins_charge = 44.27674696
    assert np.isclose(converged_results['total charge in muffin-tins'], muffin_tins_charge),\
        f"Incorrect value for total charge in muffin-tins. Expect {muffin_tins_charge}"

    total_charge = 47.0
    assert np.isclose(converged_results['total charge'], total_charge),\
        f"Incorrect value for total charge. Expect {total_charge}"

def test_loss_fxc_RPA_LF(loss_fxcRPA: np.ndarray):
    """Test results of loss function using the random-phase approximation (RPA)
     for the exchange-correlation kernel when local-field effects are included.
    """
    assert np.allclose(loss_fxcRPA[:, 0], ref_energies),\
        f"Energy grid not equivalent to reference calculations using the random-phase approximation (RPA) for the exchange-correlation kernel when local-field effects are included"
    assert np.allclose(loss_fxcRPA[:, 1], ref_loss_function_RPA_LF),\
        f"Loss function not equivalent to reference calculations using the random-phase approximation (RPA) for the exchange-correlation kernel when local-field effects are included"
    
def test_loss_fxc_RPA_NLF(loss_fxcRPA_nlf: np.ndarray):
    """Test results of loss function using the random-phase approximation (RPA)
    for the exchange-correlation kernel when local-field effects are neglected.
    """
    assert np.allclose(loss_fxcRPA_nlf[:, 0], ref_energies),\
        f"Energy grid not equivalent to reference calculations using the random-phase approximation (RPA) for the exchange-correlation kernel when local-field effects are neglected"
    assert np.allclose(loss_fxcRPA_nlf[:, 1], ref_loss_function_RPA_NLF),\
        f"Loss function not equivalent to reference calculations using the random-phase approximation (RPA) for the exchange-correlation kernel when local-field effects are neglected"
    
def test_loss_fxc_ALDA_LF(loss_fxcALDA: np.ndarray):
    """Test results of loss function using the adiabatic LDA (ALDA)
     for the exchange-correlation kernel when local-field effects are included.
    """
    assert np.allclose(loss_fxcALDA[:, 0], ref_energies),\
        f"Energy grid not equivalent to reference calculations using the adiabatic LDA (ALDA) for the exchange-correlation kernel when local-field effects are included"
    assert np.allclose(loss_fxcALDA[:, 1], ref_loss_function_ALDA_LF),\
        f"Loss function not equivalent to reference calculations using the adiabatic LDA (ALDA) for the exchange-correlation kernel when local-field effects are included"

def test_loss_fxc_ALDA_NLF(loss_fxcALDA_nlf: np.ndarray):
    """Test results of loss function using the adiabatic LDA (ALDA)
     for the exchange-correlation kernel when local-field effects are neglected.
    """
    assert np.allclose(loss_fxcALDA_nlf[:, 0], ref_energies), \
        f"Energy grid not equivalent to reference calculations using the adiabatic LDA (ALDA) for the exchange-correlation kernel when local-field effects are neglected"
    assert np.allclose(loss_fxcALDA_nlf[:, 1], ref_loss_function_ALDA_NLF),\
        f"Loss function not equivalent to reference calculations using the adiabatic LDA (ALDA) for the exchange-correlation kernel when local-field effects are neglected"

def test_dielectric_fxc_RPA_LF(dielectric_fxcRPA: np.ndarray):
    """Test results of dielectric function using the random-phase approximation (RPA)
    for the exchange-correlation kernel when local-field effects are included.
    """
    assert np.allclose(dielectric_fxcRPA[:, 0], ref_frequency),\
        f"Frequency grid not equivalent to reference calculations using the random-phase approximation (RPA) for the exchange-correlation kernel when local-field effects are included"
    assert np.allclose(dielectric_fxcRPA[:, 1], ref_real_dielectric_function_RPA_LF),\
        f"Real part of dielectric function not equivalent to reference calculations using the random-phase approximation (RPA) for the exchange-correlation kernel when local-field effects are included"
    assert np.allclose(dielectric_fxcRPA[:, 2], ref_imag_dielectric_function_RPA_LF),\
        f"Imaginary part of dielectric function not equivalent to reference calculations using the random-phase approximation (RPA) for the exchange-correlation kernel when local-field effects are included"
    assert np.allclose(dielectric_fxcRPA[:, 3], ref_real_KKT_dielectric_function_RPA_LF),\
        f"Real part of dielectric function (by Kramers-Kronig) not equivalent to reference calculations using the random-phase approximation (RPA) for the exchange-correlation kernel when local-field effects are included"

def test_dielectric_fxc_RPA_NLF(dielectric_fxcRPA_nlf: np.ndarray):
    """Test results of dielectric function using the random-phase approximation (RPA)
    for the exchange-correlation kernel when local-field effects are neglected.
    """
    assert np.allclose(dielectric_fxcRPA_nlf[:, 0], ref_frequency),\
        f"Frequency grid not equivalent to reference calculations using the random-phase approximation (RPA) for the exchange-correlation kernel when local-field effects are neglected"
    assert np.allclose(dielectric_fxcRPA_nlf[:, 1], ref_real_dielectric_function_RPA_NLF),\
        f"Real part of dielectric function not equivalent to reference calculations using the random-phase approximation (RPA) for the exchange-correlation kernel when local-field effects are neglected"
    assert np.allclose(dielectric_fxcRPA_nlf[:, 2], ref_imag_dielectric_function_RPA_NLF),\
        f"Imaginary part of dielectric function not equivalent to reference calculations using the random-phase approximation (RPA) for the exchange-correlation kernel when local-field effects are neglected"
    assert np.allclose(dielectric_fxcRPA_nlf[:, 3], ref_real_KKT_dielectric_function_RPA_NLF),\
        f"Real part of dielectric function (by Kramers-Kronig) not equivalent to reference calculations using the random-phase approximation (RPA) for the exchange-correlation kernel when local-field effects are neglected"

def test_dielectric_fxc_ALDA_LF(dielectric_fxcALDA: np.ndarray):
    """Test results of dielectric function using the adiabatic LDA (ALDA)
    for the exchange-correlation kernel when local-field effects are included.
    """
    assert np.allclose(dielectric_fxcALDA[:, 0], ref_frequency),\
        f"Frequency grid not equivalent to reference calculations using the adiabatic LDA (ALDA) for the exchange-correlation kernel when local-field effects are included"
    assert np.allclose(dielectric_fxcALDA[:, 1], ref_real_dielectric_function_ALDA_LF),\
        f"Real part of dielectric function not equivalent to reference calculations using the adiabatic LDA (ALDA) for the exchange-correlation kernel when local-field effects are included"
    assert np.allclose(dielectric_fxcALDA[:, 2], ref_imag_dielectric_function_ALDA_LF),\
        f"Imaginary part of dielectric function not equivalent to reference calculations using the adiabatic LDA (ALDA) for the exchange-correlation kernel when local-field effects are included"
    assert np.allclose(dielectric_fxcALDA[:, 3], ref_real_KKT_dielectric_function_ALDA_LF),\
        f"Real part of dielectric function (by Kramers-Kronig) not equivalent to reference calculations using then adiabatic LDA (ALDA)for the exchange-correlation kernel when local-field effects are included"

def test_dielectric_fxc_ALDA_NLF(dielectric_fxcALDA_nlf: np.ndarray):
    """Test results of dielectric function using the adiabatic LDA (ALDA)
    for the exchange-correlation kernel when local-field effects are neglected.
    """
    assert np.allclose(dielectric_fxcALDA_nlf[:, 0], ref_frequency),\
        f"Frequency grid not equivalent to reference calculations using the adiabatic LDA (ALDA) for the exchange-correlation kernel when local-field effects are neglected"
    assert np.allclose(dielectric_fxcALDA_nlf[:, 1], ref_real_dielectric_function_ALDA_NLF),\
        f"Real part of dielectric function not equivalent to reference calculations using the adiabatic LDA (ALDA) for the exchange-correlation kernel when local-field effects are neglected"
    assert np.allclose(dielectric_fxcALDA_nlf[:, 2], ref_imag_dielectric_function_ALDA_NLF),\
        f"Imaginary part of dielectric function not equivalent to reference calculations using the adiabatic LDA (ALDA) for the exchange-correlation kernel when local-field effects are neglected"
    assert np.allclose(dielectric_fxcALDA_nlf[:, 3], ref_real_KKT_dielectric_function_ALDA_NLF),\
        f"Real part of dielectric function (by Kramers-Kronig) not equivalent to reference calculations using the adiabatic LDA (ALDA)for the exchange-correlation kernel when local-field effects are neglected"