from .execute import execute
from ..infrastructure import flatten_directory, remove_ignored_files


def run_single_reference(job_dir: str, executable: str, settings):
    """
    Reference run for single test case.

    To facilitate the test suite, exciting outputs that are written to subdirectories
    are extracted and placed in the job directory.

    exciting outputs that are not used by the test suite are deleted (else one can
    be left with large binary files, etc, that are erroneously added to the repo).

    :param str job_dir: Job directory
    :param executable: exciting executable
    :param namedtuple settings: Default settings for max calculation timing, output
     file name and exciting outputs to be ignored (hence deleted) by the test suite.
    """

    success, err_mess, timing = execute(job_dir, executable, settings.main_output, settings.max_time)

    if success:
        print('%s: Run succeeded' % job_dir)
        print('Time (s): %.1f' % timing)
        flatten_directory(job_dir)
        remove_ignored_files(job_dir, settings.ignored_output)

    else:
        print('%s: Run failed' % job_dir)
        print(*(err for err in err_mess), sep='\n')
        print('Time (s): %.1f' % timing)
