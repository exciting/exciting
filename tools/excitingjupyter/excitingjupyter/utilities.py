import os
import pathlib

def get_exciting_root() -> str:
    """Get exciting's root directory.

    :return exciting root directory.
    """
    if not 'EXCITINGROOT' in os.environ.keys():
        print(f"wd in get_exciting_root: {os.getcwd()}")
        excitingroot = str(pathlib.Path(os.getcwd()).parents[3]) 
        os.environ['EXCITINGROOT'] = excitingroot
    else:
        excitingroot = os.environ['EXCITINGROOT']
    assert "exciting_smp" in os.listdir(os.path.join(excitingroot, "bin")), f"Exciting binary 'exciting_smp' was not found at exciting root path: {excitingroot}"
    return excitingroot
