import os.path
import numpy as np
from excitingtools import parser_chooser


def test_groundstate_BN(converged_results_BN: dict):
    """Test results (for cubic boron-nitride) in main output file INFO.OUT
    containing the essential information on the material system, parameters of the calculation and results
    """
   
    total_energy = -79.38465862
    assert np.isclose(converged_results_BN['Total energy'], total_energy),\
        f"Incorrect value for total energy. Expect {total_energy}"

    fermi_energy = 0.42726809
    assert np.isclose(converged_results_BN['Fermi energy'], fermi_energy),\
        f"Incorrect value for fermi energy. Expect {fermi_energy}"

    kinetic_energy = 79.03180996
    assert np.isclose(converged_results_BN['Kinetic energy'], kinetic_energy),\
        f"Incorrect value for kinetic energy. Expect {kinetic_energy}"

    coulomb_energy = -147.71206057
    assert np.isclose(converged_results_BN['Coulomb energy'], coulomb_energy),\
        f"Incorrect value for coulomb energy. Expect {coulomb_energy}"

    exchange_energy = -10.2100207
    assert np.isclose(converged_results_BN['Exchange energy'], exchange_energy),\
        f"Incorrect value for exchange energy. Expect {exchange_energy}"

    correlation_energy = -0.49438731
    assert np.isclose(converged_results_BN['Correlation energy'], correlation_energy),\
        f"Incorrect value for correlation energy. Expect {correlation_energy}"

    fermi_energy_dos = 0.0
    assert np.isclose(converged_results_BN['DOS at Fermi energy (states/Ha/cell)'], fermi_energy_dos),\
        f"Incorrect value for DOS at Fermi energy (states/Ha/cell). Expect {fermi_energy_dos}"

    charge_core_electrons = 4.0
    assert np.isclose(converged_results_BN['core'], charge_core_electrons),\
        f"Incorrect value for charge of core electrons. Expect {charge_core_electrons}"

    core_leakage_charge = 0.0014585
    assert np.isclose(converged_results_BN['core leakage'], core_leakage_charge),\
        f"Incorrect value for core leakage charge. Expect {core_leakage_charge}"

    charge_valance_electrons = 8.0
    assert np.isclose(converged_results_BN['valence'], charge_valance_electrons),\
        f"Incorrect value for charge of valence electrons. Expect {charge_valance_electrons}"

    interstitial_region_charge = 3.23877764
    assert np.isclose(converged_results_BN['interstitial'], interstitial_region_charge),\
        f"Incorrect value for charge in interstitial region. Expect {interstitial_region_charge}"

    muffin_tins_charge = 8.76122236
    assert np.isclose(converged_results_BN['total charge in muffin-tins'], muffin_tins_charge),\
        f"Incorrect value for total charge in muffin-tins. Expect {muffin_tins_charge}"

    total_charge = 12.0
    assert np.isclose(converged_results_BN['total charge'], total_charge),\
        f"Incorrect value for total charge. Expect {total_charge}"
    
    estimated_gap = 0.15663379
    assert np.isclose(converged_results_BN['Estimated fundamental gap'], estimated_gap),\
        f"Incorrect value for Estimated fundamental gap. Expect {estimated_gap}"
    
    
def test_dielectric_BN(epsilon_results_BN: dict, reference: str):
    """Test results of dielectric function calculations for cubic boron-nitride
    """
    
    epsilon_reference_BN = parser_chooser(os.path.join(reference))

    assert np.allclose(epsilon_results_BN['frequency'], epsilon_reference_BN['frequency']), \
        f"Frequency grid not equivalent to reference calculation for dielectric function calculations for cubic boron-nitride"

    assert np.allclose(epsilon_results_BN['real_oscillator_strength'], epsilon_reference_BN['real_oscillator_strength']), \
        f"Real part of dielectric function not equivalent to reference calculation for cubic boron-nitride"

    assert np.allclose(epsilon_results_BN['imag_oscillator_strength'], epsilon_reference_BN['imag_oscillator_strength']), \
        f"Imaginary part of dielectric function not equivalent to reference calculation for cubic boron-nitride"
                        
    assert np.allclose(epsilon_results_BN['real_oscillator_strength_kkt'], epsilon_reference_BN['real_oscillator_strength_kkt']),\
         f"Real part of dielectric function (by Kramers-Kronig) not equivalent to reference calculation for cubic boron-nitride"

    
def test_loss_BN(loss_results_BN: dict, reference: str):
    """Test results of loss function calculations for cubic boron-nitride
    """
    
    loss_reference_BN = parser_chooser(os.path.join(reference))

    assert np.allclose(loss_results_BN['frequency'], loss_reference_BN['frequency']), \
        f"Frequency grid not equivalent to reference calculation for loss function calculations for cubic boron-nitride"

    assert np.allclose(loss_results_BN['real_oscillator_strength'], loss_reference_BN['real_oscillator_strength']), \
        f"Real part of loss function not equivalent to reference calculation for cubic boron-nitride"

    assert np.allclose(loss_results_BN['imag_oscillator_strength'], loss_reference_BN['imag_oscillator_strength']), \
        f"Imaginary part of loss function not equivalent to reference calculation for cubic boron-nitride"     
    

def test_groundstate_TiO2(converged_results_TiO2: dict):
    """Test results (for rutile titanium dioxide) in main output file INFO.OUT
    containing the essential information on the material system, parameters of the calculation and results
    """
    
    total_energy = -2006.17125266
    assert np.isclose(converged_results_TiO2['Total energy'], total_energy),\
        f"Incorrect value for total energy. Expect {total_energy}"

    fermi_energy = 0.21453506
    assert np.isclose(converged_results_TiO2['Fermi energy'], fermi_energy),\
        f"Incorrect value for fermi energy. Expect {fermi_energy}"

    kinetic_energy = 2021.14376199
    assert np.isclose(converged_results_TiO2['Kinetic energy'], kinetic_energy),\
        f"Incorrect value for kinetic energy. Expect {kinetic_energy}"

    coulomb_energy = -3911.79183723
    assert np.isclose(converged_results_TiO2['Coulomb energy'], coulomb_energy),\
        f"Incorrect value for coulomb energy. Expect {coulomb_energy}"

    exchange_energy = -112.03793785
    assert np.isclose(converged_results_TiO2['Exchange energy'], exchange_energy),\
        f"Incorrect value for exchange energy. Expect {exchange_energy}"

    correlation_energy = -3.48523958
    assert np.isclose(converged_results_TiO2['Correlation energy'], correlation_energy),\
        f"Incorrect value for correlation energy. Expect {correlation_energy}"

    fermi_energy_dos = 0.0
    assert np.isclose(converged_results_TiO2['DOS at Fermi energy (states/Ha/cell)'], fermi_energy_dos),\
        f"Incorrect value for DOS at Fermi energy (states/Ha/cell). Expect {fermi_energy_dos}"

    charge_core_electrons = 28.0
    assert np.isclose(converged_results_TiO2['core'], charge_core_electrons),\
        f"Incorrect value for charge of core electrons. Expect {charge_core_electrons}"

    core_leakage_charge = 1.11e-06
    assert np.isclose(converged_results_TiO2['core leakage'], core_leakage_charge),\
        f"Incorrect value for core leakage charge. Expect {core_leakage_charge}"

    charge_valance_electrons = 48.0
    assert np.isclose(converged_results_TiO2['valence'], charge_valance_electrons),\
        f"Incorrect value for charge of valence electrons. Expect {charge_valance_electrons}"

    interstitial_region_charge = 7.2840741
    assert np.isclose(converged_results_TiO2['interstitial'], interstitial_region_charge),\
        f"Incorrect value for charge in interstitial region. Expect {interstitial_region_charge}"

    muffin_tins_charge = 68.71592584
    assert np.isclose(converged_results_TiO2['total charge in muffin-tins'], muffin_tins_charge),\
        f"Incorrect value for total charge in muffin-tins. Expect {muffin_tins_charge}"

    total_charge = 76.0
    assert np.isclose(converged_results_TiO2['total charge'], total_charge),\
        f"Incorrect value for total charge. Expect {total_charge}"
    
    estimated_gap = 0.06638281
    assert np.isclose(converged_results_TiO2['Estimated fundamental gap'], estimated_gap),\
        f"Incorrect value for Estimated fundamental gap. Expect {estimated_gap}"
    
    
def test_dielectric_TiO2(epsilon_results_TiO2: dict, reference: str):
    """Test results of dielectric function calculations for rutile titanium dioxide
    """
    
    epsilon_reference_TiO2 = parser_chooser(os.path.join(reference))

    assert np.allclose(epsilon_results_TiO2['frequency'], epsilon_reference_TiO2['frequency']), \
        f"Frequency grid not equivalent to reference calculation for dielectric function calculations for rutile titanium dioxide"

    assert np.allclose(epsilon_results_TiO2['real_oscillator_strength'], epsilon_reference_TiO2['real_oscillator_strength']), \
        f"Real part of dielectric function not equivalent to reference calculation for rutile titanium dioxide"

    assert np.allclose(epsilon_results_TiO2['imag_oscillator_strength'], epsilon_reference_TiO2['imag_oscillator_strength']), \
        f"Imaginary part of dielectric function not equivalent to reference calculation for rutile titanium dioxide"
                        
    assert np.allclose(epsilon_results_TiO2['real_oscillator_strength_kkt'], epsilon_reference_TiO2['real_oscillator_strength_kkt']), \
        f"Real part of dielectric function (by Kramers-Kronig) not equivalent to reference calculation for rutile titanium dioxide"

    
def test_loss_TiO2(loss_results_TiO2: dict, reference: str):
    """Test results of loss function calculations for rutile titanium dioxide
    """
    
    loss_reference_TiO2 = parser_chooser(os.path.join(reference))

    assert np.allclose(loss_results_TiO2['frequency'], loss_reference_TiO2['frequency']), \
        f"Frequency grid not equivalent to reference calculation for loss function calculations for rutile titanium dioxide"

    assert np.allclose(loss_results_TiO2['real_oscillator_strength'], loss_reference_TiO2['real_oscillator_strength']), \
        f"Real part of loss function not equivalent to reference calculation for rutile titanium dioxide"

    assert np.allclose(loss_results_TiO2['imag_oscillator_strength'], loss_reference_TiO2['imag_oscillator_strength']), \
        f"Imaginary part of loss function not equivalent to reference calculation for rutile titanium dioxide"     