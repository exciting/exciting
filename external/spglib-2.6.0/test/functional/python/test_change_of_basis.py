from __future__ import annotations

import numpy as np
from spglib import get_symmetry_from_database


def test_change_of_basis(crystal_data_dataset):
    crystal_data = crystal_data_dataset["crystal_data"]
    dataset = crystal_data_dataset["dataset"]
    symprec = crystal_data_dataset["symprec"]
    std_pos = dataset.std_positions
    tmat = dataset.transformation_matrix
    orig_shift = dataset.origin_shift
    lat = np.dot(crystal_data.cell[0].T, np.linalg.inv(tmat))
    pos = np.dot(crystal_data.cell[1], tmat.T) + orig_shift
    for p in pos:
        diff = std_pos - p
        diff -= np.rint(diff)
        diff = np.dot(diff, lat.T)
        delta = np.sqrt((diff**2).sum(axis=1))
        indices = np.where(delta < symprec)[0]
        assert len(indices) == 1


def test_std_symmetry(crystal_data_dataset):
    dataset = crystal_data_dataset["dataset"]
    symmetry = get_symmetry_from_database(dataset.hall_number)
    std_pos = dataset.std_positions

    # for r, t in zip(symmetry['rotations'], symmetry['translations']):
    #     for rp in (np.dot(std_pos, r.T) + t):
    #         diff = std_pos - rp
    #         diff -= np.rint(diff)
    #         num_match = len(np.where(abs(diff).sum(axis=1) < 1e-3)[0])
    #         self.assertEqual(num_match, 1, msg="%s" % fname)

    # Equivalent above by numpy hack
    # 15 sec on macOS 2.3 GHz Intel Core i5 (4times faster than above)
    rot = symmetry["rotations"]
    trans = symmetry["translations"]
    # (n_sym, 3, n_atom)
    rot_pos = np.dot(rot, std_pos.T) + trans[:, :, None]
    for p in std_pos:
        diff = rot_pos - p[None, :, None]
        diff -= np.rint(diff)
        num_match = (abs(diff).sum(axis=1) < 1e-3).sum(axis=1)
        assert all(num_match == 1)


def test_std_rotation(crystal_data_dataset):
    crystal_data = crystal_data_dataset["crystal_data"]
    dataset = crystal_data_dataset["dataset"]
    symprec = crystal_data_dataset["symprec"]
    std_lat = dataset.std_lattice
    tmat = dataset.transformation_matrix
    lat = np.dot(crystal_data.cell[0].T, np.linalg.inv(tmat))
    lat_rot = np.dot(dataset.std_rotation_matrix, lat)
    np.testing.assert_allclose(
        std_lat,
        lat_rot.T,
        atol=symprec,
    )
