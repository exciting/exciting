from ..runner.execute import execute
from ..io.file_system import flatten_directory, remove_ignored_files


def run_single_reference(job_dir: str, execution_cmd: str, settings, my_env):
    """
    Reference run for single test case.

    To facilitate the test suite, exciting outputs that are written to subdirectories
    are extracted and placed in the job directory.

    exciting outputs that are not used by the test suite are deleted (else one can
    be left with large binary files, etc, that are erroneously added to the repo).

    :param str job_dir: Job directory.
    :param execution_cmd: Execution command.
    :param namedtuple settings: Default settings for max calculation timing, output
     file name and exciting outputs to be ignored (hence deleted) by the test suite.
    :param my_env: An instance of the shell environment.
    """
    success, err_mess, timing = execute(job_dir, execution_cmd, settings.max_time, my_env=my_env)

    if success:
        print('%s: Run exited cleanly' % job_dir)
        print('Time (s): %.1f' % timing)
        flatten_directory(job_dir)
        remove_ignored_files(job_dir, settings.ignored_output)

    else:
        print('%s: Run failed' % job_dir)
        print(*(err for err in err_mess), sep='\n')
        print('Time (s): %.1f' % timing)
