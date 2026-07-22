from pathlib import Path

import pytest
from excitingscripts.setup.excited_state_species_files import (
    DATA_DIR,
    DEFAULT_LO_PATH,
    DEFAULT_OUTPUT_PATH,
    DEFAULT_XML_PATH,
    determine_required_states,
    main,
    optimize_species_for_excited_states,
    parse_lo_recommendation,
    validate_species_file,
)
from excitingtools.species.species_file import SpeciesFile
from excitingtools.utils.test_utils import MockFile


@pytest.fixture
def ground_state_species_mock(tmp_path) -> MockFile:
    """ Mock ground-state 'Si.xml' species file (l=0 and l=1 only).
    """
    species_str = """<?xml version="1.0" ?>
<spdb>
	<sp chemicalSymbol="Si" name="silicon" z="-14.0" mass="51196.73454">
		<muffinTin rmin="1e-06" radius="2.0" rinf="24.976" radialmeshPoints="600"/>
		<atomicState n="1" l="0" kappa="1" occ="2.0" core="true"/>
		<atomicState n="2" l="0" kappa="1" occ="2.0" core="false"/>
		<atomicState n="2" l="1" kappa="1" occ="2.0" core="false"/>
		<atomicState n="2" l="1" kappa="2" occ="4.0" core="false"/>
		<atomicState n="3" l="0" kappa="1" occ="2.0" core="false"/>
		<atomicState n="3" l="1" kappa="1" occ="1.0" core="false"/>
		<atomicState n="3" l="1" kappa="2" occ="1.0" core="false"/>
		<basis>
			<default type="lapw" trialEnergy="0.15" searchE="false"/>
			<custom l="0" type="lapw" n="3" searchE="false"/>
			<custom l="1" type="lapw" n="3" searchE="false"/>
			<lo l="0">
				<wf matchingOrder="0" searchE="false" n="3"/>
				<wf matchingOrder="0" searchE="false" n="2"/>
			</lo>
			<lo l="1">
				<wf matchingOrder="0" searchE="false" n="3"/>
				<wf matchingOrder="0" searchE="false" n="2"/>
			</lo>
			<lo l="0">
				<wf matchingOrder="0" searchE="false" n="2"/>
				<wf matchingOrder="1" searchE="false" n="2"/>
			</lo>
			<lo l="0">
				<wf matchingOrder="0" searchE="false" n="3"/>
				<wf matchingOrder="1" searchE="false" n="3"/>
			</lo>
			<lo l="1">
				<wf matchingOrder="0" searchE="false" n="2"/>
				<wf matchingOrder="1" searchE="false" n="2"/>
			</lo>
			<lo l="1">
				<wf matchingOrder="0" searchE="false" n="3"/>
				<wf matchingOrder="1" searchE="false" n="3"/>
			</lo>
		</basis>
	</sp>
</spdb>
"""
    species_file = tmp_path / "Si.xml"
    species_file.write_text(species_str)
    return MockFile(species_file, species_str)


@pytest.fixture
def lo_recommendation_mock(tmp_path) -> MockFile:
    """ Mock 'LO_RECOMMENDATION.OUT' with recommendations for l=0..8, spanning enough 'n' to
    exercise the full energy_threshold=80.0 Ha cutoff for every channel.
    """
    lo_recommendation_str = """# Recommended linearization energies computet with Wigner-Seitz rules.
 --------------------------------------------------------------------
 #  n_species:  1
 # n_l-channels:  9
 # n_nodes: 11

 # species: Si, l :  0
 # nodes   n        trial energy
       0   1    -65.093766647505
       1   2     -4.688526320986
       2   3      0.485021234065
       3   4      5.211392739027
       4   5     12.977737805228
       5   6     23.346022327098
       6   7     36.189576708965
       7   8     51.442359263906
       8   9     69.061394085392
       9  10     89.017161125832
      10  11    111.288623673197

 # species: Si, l :  1
 # nodes   n        trial energy
       0   2     -3.101736075221
       1   3      0.921064941142
       2   4      5.413422040881
       3   5     12.644836187708
       4   6     22.315631543845
       5   7     34.350487259126
       6   8     48.715880164301
       7   9     65.391188013595
       8  10     84.362569808554
       9  11    105.619976898581
      10  12    129.155783783329

 # species: Si, l :  2
 # nodes   n        trial energy
       0   3      1.387055391814
       1   4      5.041019390005
       2   5     11.190751677305
       3   6     19.739310489307
       4   7     30.647658191658
       5   8     43.889515765000
       6   9     59.444189224985
       7  10     77.296671715763
       8  11     97.436040600560
       9  12    119.854223058443
      10  13    144.545099296410

 # species: Si, l :  3
 # nodes   n        trial energy
       0   4      3.162527632538
       1   5      8.450768720686
       2   6     16.009727893100
       3   7     25.880320708435
       4   8     38.050535206396
       5   9     52.515748453878
       6  10     69.269889767019
       7  11     88.306499466767
       8  12    109.619493647007
       9  13    133.203584847771
      10  14    159.054341506809

 # species: Si, l :  4
 # nodes   n        trial energy
       0   5      5.014761856696
       1   6     11.754337727125
       2   7     20.636437158384
       3   8     31.782792581634
       4   9     45.197844019409
       5  10     60.882003964334
       6  11     78.835104395074
       7  12     99.055830572911
       8  13    121.542087045472
       9  14    146.291431908831
      10  15    173.301390871190

 # species: Si, l :  5
 # nodes   n        trial energy
       0   6      7.072717988966
       1   7     15.204929271203
       2   8     25.363486786148
       3   9     37.747956773595
       4  10     52.381435013247
       5  11     69.268310556532
       6  12     88.410745731329
       7  13    109.809667599090
       8  14    133.465125396444
       9  15    159.376542715683
      10  16    187.542960845926

 # species: Si, l :  6
 # nodes   n        trial energy
       0   7      9.364438623209
       1   8     18.873220724371
       2   9     30.287887250145
       3  10     43.890746579252
       4  11     59.724898483463
       5  12     77.800770260432
       6  13     98.122686564240
       7  14    120.692935461160
       8  15    145.512749460499
       9  16    172.582656521678
      10  17    201.902694602855

 # species: Si, l :  7
 # nodes   n        trial energy
       0   8     11.899151372175
       1   9     22.785866459404
       2  10     35.449241555894
       3  11     50.261487833095
       4  12     67.286367051900
       5  13     86.542123308346
       6  14    108.036008199369
       7  15    131.771742290281
       8  16    157.751489351977
       9  17    185.976536291223
      10  18    216.447597213935

 # species: Si, l :  8
 # nodes   n        trial energy
       0   9     14.680299694182
       1  10     26.954517246801
       2  11     40.866068270314
       3  12     56.884734845433
       4  13     75.095643034762
       5  14     95.526032226256
       6  15    118.186933655405
       7  16    143.083881895688
       8  17    170.220062053112
       9  18    199.597455169570
      10  19    231.217323093414
"""
    lo_file = tmp_path / "LO_RECOMMENDATION.OUT"
    lo_file.write_text(lo_recommendation_str)
    return MockFile(lo_file, lo_recommendation_str)


def _basis_signature(species_obj: SpeciesFile):
    """ Order-independent signature of a species file's basis: custom LAPWs (as (l, n) pairs) and
    local orbitals (as (l, sorted n's, sorted matchingOrders) triplets).
    """
    custom = {(c["l"], c["n"]) for c in species_obj.basis["custom"]}
    lo = {
        (lo["l"], tuple(sorted(wf["n"] for wf in lo["wf"])), tuple(sorted(wf["matchingOrder"] for wf in lo["wf"])))
        for lo in species_obj.basis["lo"]
    }
    return custom, lo


def test_optimize_species_for_excited_states(ground_state_species_mock, lo_recommendation_mock, tmp_path):
    output_dir = tmp_path / "excited_state_species"

    optimize_species_for_excited_states(
        species_name="Si",
        energy_threshold=80.0,
        path_xml=tmp_path,
        path_lo=tmp_path,
        path_out=output_dir,
        max_matching_order=1,
    )

    output_file = output_dir / "Si_excited.xml"
    assert output_file.exists()

    result_obj = SpeciesFile.from_file(output_file)
    custom_signature, lo_signature = _basis_signature(result_obj)

    # New l-channels 2..8 must get a custom LAPW at n = l + 1, in addition to the
    # ground-state custom LAPWs already present at l=0, l=1.
    expected_custom = {(0, 3), (1, 3), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (8, 9)}
    assert custom_signature == expected_custom, "Custom APWs are incorrect."

    # Highest 'n' reached per l-channel, given energy_threshold=80.0 Ha (derived from the
    # trial energies above: the last state <= 80 Ha, and the first state > 80 Ha, per l).
    highest_n = {0: 9, 1: 9, 2: 10, 3: 10, 4: 11, 5: 11, 6: 12, 7: 12, 8: 13}
    lowest_n = {0: 4, 1: 4, 2: 3, 3: 4, 4: 5, 5: 6, 6: 7, 7: 8, 8: 9}

    for l, n_max in highest_n.items():
        n_min = lowest_n[l]
        # (0, 1) high-energy LO must exist for every required n in this channel.
        for n in range(n_min, n_max + 1):
            assert (l, (n, n), (0, 1)) in lo_signature, f"missing (0,1) LO for l={l}, n={n}"
        # Bridging (0, 0) LOs must connect every adjacent pair, except the very first state of a
        # brand-new l-channel (l=2..8), which has no lower state to bridge from.
        bridge_start = n_min + 1 if l >= 2 else n_min
        for n in range(bridge_start, n_max + 1):
            assert (l, (n - 1, n), (0, 0)) in lo_signature, f"missing bridging LO for l={l}, n=({n - 1},{n})"
        # No LOs above n_max for this channel (i.e. the threshold cutoff was respected).
        assert not any(sig[0] == l and n_max + 1 in sig[1] for sig in lo_signature)

    # No matching order beyond 1 should have been added anywhere (max_matching_order=1).
    assert all(max(sig[2]) <= 1 for sig in lo_signature), (
        "No matching order beyond 1 should have been added anywhere (max_matching_order=1)"
    )


def test_default_paths_are_in_data_directory():
    assert DEFAULT_XML_PATH == DATA_DIR / "ground_state_species"
    assert DEFAULT_LO_PATH == DATA_DIR / "lo_recommendations"
    assert DEFAULT_OUTPUT_PATH == DATA_DIR / "excited_state_species"


def test_parse_lo_recommendation_selects_requested_species(tmp_path):
    recommendation_file = tmp_path / "LO_RECOMMENDATION.OUT"
    recommendation_file.write_text(
        """# species: Si, l :  0
0 1 -1.0
1 2 2.5
# species: C, l :  0
0 1 -2.0
"""
    )

    assert parse_lo_recommendation(recommendation_file, "Si") == {
        0: [(0, 1, -1.0), (1, 2, 2.5)]
    }


def test_energy_cutoff_selects_only_states_up_to_threshold(ground_state_species_mock):
    species_obj = SpeciesFile.from_file(ground_state_species_mock.file)
    recommendations = {0: [(3, 4, 5.0), (4, 5, 10.0), (5, 6, 20.0)]}

    assert determine_required_states(species_obj, recommendations, 10.0) == {0: 2}


@pytest.mark.parametrize("max_matching_order", [-1, 4])
def test_invalid_max_matching_order_is_rejected(
    max_matching_order, ground_state_species_mock, lo_recommendation_mock, tmp_path
):
    with pytest.raises(ValueError, match="max_matching_order must be between 0 and 3"):
        optimize_species_for_excited_states(
            species_name="Si",
            energy_threshold=80.0,
            path_xml=tmp_path,
            path_lo=tmp_path,
            path_out=tmp_path / "out",
            max_matching_order=max_matching_order,
        )


def test_cli_rejects_invalid_max_matching_order(monkeypatch):
    monkeypatch.setattr(
        "sys.argv",
        ["excited_state_species_files", "Si", "-e", "80", "--max-matching-order", "4"],
    )

    with pytest.raises(SystemExit, match="2"):
        main()


@pytest.mark.parametrize("basis_kind", ["custom", "lo"])
def test_trial_energy_species_files_are_rejected(ground_state_species_mock, basis_kind):
    species_obj = SpeciesFile.from_file(ground_state_species_mock.file)
    if basis_kind == "custom":
        species_obj.basis["custom"][0]["trialEnergy"] = 1.0
    else:
        species_obj.basis["lo"][0]["wf"][0]["trialEnergy"] = 1.0

    with pytest.raises(ValueError, match="explicit trialEnergy"):
        validate_species_file(species_obj)


def test_existing_custom_lapw_gets_missing_first_helo(
    ground_state_species_mock, lo_recommendation_mock, tmp_path
):
    species_obj = SpeciesFile.from_file(ground_state_species_mock.file)
    species_obj.basis["custom"].append({"l": 2, "type": "lapw", "n": 3, "searchE": False})
    species_obj.write(ground_state_species_mock.file)

    output_dir = tmp_path / "excited_state_species"
    optimize_species_for_excited_states(
        species_name="Si",
        energy_threshold=80.0,
        path_xml=tmp_path,
        path_lo=tmp_path,
        path_out=output_dir,
        max_matching_order=1,
    )

    result_obj = SpeciesFile.from_file(output_dir / "Si_excited.xml")
    custom_signature, lo_signature = _basis_signature(result_obj)

    assert (2, 3) in custom_signature
    assert (2, (2, 3), (0, 0)) not in lo_signature
    assert (2, (3, 3), (0, 1)) in lo_signature
    assert (2, (3, 4), (0, 0)) in lo_signature
    assert (2, (4, 4), (0, 1)) in lo_signature


def test_existing_custom_apw_lo_is_not_duplicated(
    ground_state_species_mock, lo_recommendation_mock, tmp_path
):
    species_obj = SpeciesFile.from_file(ground_state_species_mock.file)
    species_obj.basis["custom"].append({"l": 2, "type": "apw+lo", "n": 3, "searchE": False})
    species_obj.write(ground_state_species_mock.file)

    output_dir = tmp_path / "excited_state_species"
    optimize_species_for_excited_states(
        species_name="Si",
        energy_threshold=80.0,
        path_xml=tmp_path,
        path_lo=tmp_path,
        path_out=output_dir,
        max_matching_order=1,
    )

    result_obj = SpeciesFile.from_file(output_dir / "Si_excited.xml")
    _custom_signature, lo_signature = _basis_signature(result_obj)

    assert (2, (2, 3), (0, 0)) not in lo_signature
    assert (2, (3, 3), (0, 1)) not in lo_signature
    assert (2, (3, 4), (0, 0)) in lo_signature
    assert (2, (4, 4), (0, 1)) in lo_signature
