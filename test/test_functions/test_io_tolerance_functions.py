import pytest

from ..src.io.tolerances import list_tolerance_files_in_directory, update_wildcard_files_under_test, \
    update_wildcard_tolerances


def test_list_tolerance_files_in_directory(tmpdir):
    """
    tmpdir is a temporary directory that exists for the scope of the test.
    """
    file1 = tmpdir / "file1.txt"
    file2 = tmpdir / "junk.json"
    file3 = tmpdir / "tolerance_groundstate.json"

    # Required, else the files won't exist (would be preferable to 'touch')
    file1.write("Arbitrary text")
    file2.write("Arbitrary text")
    file3.write("Arbitrary text")

    tolerance_files = list_tolerance_files_in_directory(tmpdir)
    assert tolerance_files == ['tolerance_groundstate.json']


def test_update_wildcard_files_under_test():
    """
    Test wildcards get replaced
    """
    files_names = ['file_01.out', 'file_01-abc.out']
    files_under_test = ['file_??.out']

    files_under_test_test = update_wildcard_files_under_test(files_under_test, files_names)

    assert files_under_test_test == ['file_01.out', 'file_01-abc.out'], \
        'Wildcards not correctly replaced in `files_under_test`'


def test_update_wildcard_files_under_test_null_behaviour():
    """
    Test files not specified in `files_under_test` are not returned by `update_wildcard_files_under_test`
    """
    files_names = ['some_other_file.out', 'specified_file.out']
    files_under_test = ['output_of_interest.out', 'specified_file.out']

    files_under_test_test = update_wildcard_files_under_test(files_under_test, files_names)

    assert files_under_test_test == ['specified_file.out'], \
        'Files not specified `files_under_test` should be absent from returned list'


def test_update_wildcard_tolerances():
    """
    Test update_wildcard_tolerances
    """
    tol1 = {'a': 5.0}
    tol2 = {'c': 47.0}

    file_names = ['file_01.out', 'file_01-abc.out', 'other_5.out', 'other_42.out']
    tolerances = {'file_??.out': tol1, 'other_5.out': tol2}

    tolerances_test = update_wildcard_tolerances(tolerances, file_names)
    tolerances_ref = {'file_01.out': tol1, 'file_01-abc.out': tol1, 'other_5.out': tol2}

    assert list(tolerances_test.keys()) == list(tolerances_ref.keys()), "Wildcards not correctly replaced"
    assert tolerances_test['file_01.out'] == tol1, "file_??.out should be assigned tol1"
    assert tolerances_test['file_01-abc.out'] == tol1, "file_??.out should be assigned tol1"
    assert tolerances_test['other_5.out'] == tol2, "other_5.out tolerances should be unchanged"


def test_update_wildcard_tolerances_erroneous_file_name():
    """
    Test update_wildcard_tolerances with an erroneous file name - containing a wildcard
    """
    file_name = ['file_01.out', 'file_??.out']
    tolerances = {'file_01.out': {'a': 1.0}, 'file_??.out': {'b': 1.1}}
    with pytest.raises(ValueError) as error_info:
        tolerances_test = update_wildcard_tolerances(tolerances, file_name)
    assert error_info.value.args[0] == "Wild card string present in file name: file_??.out"
