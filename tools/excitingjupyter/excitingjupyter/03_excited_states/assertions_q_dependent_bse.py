import os.path
import numpy as np
from excitingtools import parser_chooser


def test_groundstate(converged_results: dict):
    """Test results for groundstate calculation in main output file INFO.OUT.
    """

    total_energy = -107.19580043
    assert np.isclose(converged_results['Total energy'], total_energy), \
        f"Incorrect value for total energy. Expect {total_energy}"

    fermi_energy = 0.12491495
    assert np.isclose(converged_results['Fermi energy'], fermi_energy), \
        f"Incorrect value for fermi energy. Expect {fermi_energy}"

    kinetic_energy = 106.93717566
    assert np.isclose(converged_results['Kinetic energy'], kinetic_energy), \
        f"Incorrect value for kinetic energy. Expect {kinetic_energy}"

    coulomb_energy = -202.15690891
    assert np.isclose(converged_results['Coulomb energy'], coulomb_energy), \
        f"Incorrect value for coulomb energy. Expect {coulomb_energy}"

    exchange_energy = -11.51315050
    assert np.isclose(converged_results['Exchange energy'], exchange_energy), \
        f"Incorrect value for exchange energy. Expect {exchange_energy}"

    correlation_energy = -0.46291668
    assert np.isclose(converged_results['Correlation energy'], correlation_energy), \
        f"Incorrect value for correlation energy. Expect {correlation_energy}"

    fermi_energy_dos = 0.00000000
    assert np.isclose(converged_results['DOS at Fermi energy (states/Ha/cell)'], fermi_energy_dos), \
        f"Incorrect value for DOS at Fermi energy (states/Ha/cell). Expect {fermi_energy_dos}"

    charge_core_electrons = 2.00000000
    assert np.isclose(converged_results['core'], charge_core_electrons), \
        f"Incorrect value for charge of core electrons. Expect {charge_core_electrons}"

    core_leakage_charge = 0.00000007
    assert np.isclose(converged_results['core leakage'], core_leakage_charge), \
        f"Incorrect value for core leakage charge. Expect {core_leakage_charge}"

    charge_valance_electrons = 10.00000000
    assert np.isclose(converged_results['valence'], charge_valance_electrons), \
        f"Incorrect value for charge of valence electrons. Expect {charge_valance_electrons}"

    interstitial_region_charge = 2.03037477
    assert np.isclose(converged_results['interstitial'], interstitial_region_charge), \
        f"Incorrect value for charge in interstitial region. Expect {interstitial_region_charge}"

    muffin_tins_charge = 9.96962523
    assert np.isclose(converged_results['total charge in muffin-tins'], muffin_tins_charge), \
        f"Incorrect value for total charge in muffin-tins. Expect {muffin_tins_charge}"

    total_charge = 12.00000000
    assert np.isclose(converged_results['total charge'], total_charge), \
        f"Incorrect value for total charge. Expect {total_charge}"

    estimated_gap = 0.32832158
    assert np.isclose(converged_results['Estimated fundamental gap'], estimated_gap), \
        f"Incorrect value for Estimated fundamental gap. Expect {estimated_gap}"


def test_dielectric_optical_bse(epsilon_optical_bse_results: dict, reference: str):
    """Test results of dielectric function in optical BSE calculations.
    """

    epsilon_optical_bse_reference = parser_chooser(os.path.join(reference))

    assert np.allclose(epsilon_optical_bse_results['frequency'], epsilon_optical_bse_reference['frequency']), \
        f"Frequency grid not equivalent to reference calculation in optical BSE dielectric function calculations"

    assert np.allclose(epsilon_optical_bse_results['real_oscillator_strength'],
                       epsilon_optical_bse_reference['real_oscillator_strength'], rtol = 0.0, atol = 1e-5), \
        f"Real part of dielectric function not equivalent to reference calculation in optical BSE calculations"

    assert np.allclose(epsilon_optical_bse_results['imag_oscillator_strength'],
                       epsilon_optical_bse_reference['imag_oscillator_strength'], rtol = 0.0, atol = 1e-5), \
        f"Imaginary part of dielectric function not equivalent to reference calculation in optical BSE calculations"

    assert np.allclose(epsilon_optical_bse_results['real_oscillator_strength_kkt'],
                       epsilon_optical_bse_reference['real_oscillator_strength_kkt'], rtol = 0.0, atol = 1e-5), \
        f"Real part of dielectric function (by Kramers-Kronig) not equivalent to reference calculation in optical BSE calculations"


def test_loss_optical_bse(loss_optical_bse_results: dict, reference: str):
    """Test results of loss function in optical BSE calculations.
    """

    loss_optical_bse_reference = parser_chooser(os.path.join(reference))

    assert np.allclose(loss_optical_bse_results['frequency'], loss_optical_bse_reference['frequency']), \
        f"Frequency grid not equivalent to reference calculation in optical BSE loss function calculations"

    assert np.allclose(loss_optical_bse_results['real_oscillator_strength'], loss_optical_bse_reference['real_oscillator_strength']), \
        f"Real part of loss function not equivalent to reference calculation in optical BSE calculations"

    assert np.allclose(loss_optical_bse_results['imag_oscillator_strength'], loss_optical_bse_reference['imag_oscillator_strength']), \
        f"Imaginary part of loss function not equivalent to reference calculation in optical BSE calculations"


def test_dielectric_xray_bse(epsilon_xray_bse_results: dict, reference: str):
    """Test results of dielectric function in X-ray BSE calculations.
    """

    epsilon_xray_bse_reference = parser_chooser(os.path.join(reference))

    assert np.allclose(epsilon_xray_bse_results['frequency'], epsilon_xray_bse_reference['frequency']), \
        f"Frequency grid not equivalent to reference calculation in X-ray BSE dielectric function calculations"

    assert np.allclose(epsilon_xray_bse_results['real_oscillator_strength'],
                       epsilon_xray_bse_reference['real_oscillator_strength'], rtol = 0.0, atol = 1e-5), \
        f"Real part of dielectric function not equivalent to reference calculation in X-ray BSE calculations"

    assert np.allclose(epsilon_xray_bse_results['imag_oscillator_strength'],
                       epsilon_xray_bse_reference['imag_oscillator_strength'], rtol = 0.0, atol = 1e-5), \
        f"Imaginary part of dielectric function not equivalent to reference calculation in X-ray BSE calculations"

    assert np.allclose(epsilon_xray_bse_results['real_oscillator_strength_kkt'],
                       epsilon_xray_bse_reference['real_oscillator_strength_kkt'], rtol = 0.0, atol = 1e-5), \
        f"Real part of dielectric function (by Kramers-Kronig) not equivalent to reference calculation in X-ray BSE calculations"


def test_loss_xray_bse(loss_xray_bse_results: dict, reference: str):
    """Test results of loss function in X-ray BSE calculations.
    """

    loss_xray_bse_reference = parser_chooser(os.path.join(reference))

    assert np.allclose(loss_xray_bse_results['frequency'], loss_xray_bse_reference['frequency']), \
        f"Frequency grid not equivalent to reference calculation in X-ray BSE loss function calculations"

    assert np.allclose(loss_xray_bse_results['real_oscillator_strength'], loss_xray_bse_reference['real_oscillator_strength']), \
        f"Real part of loss function not equivalent to reference calculation in X-ray BSE calculations"

    assert np.allclose(loss_xray_bse_results['imag_oscillator_strength'], loss_xray_bse_reference['imag_oscillator_strength']), \
        f"Imaginary part of loss function not equivalent to reference calculation in X-ray BSE calculations"

