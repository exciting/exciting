"""Wrapper for command-line grep."""

import subprocess
import warnings
from pathlib import Path
from typing import Optional, Union

path_type = Union[Path, str]


def grep(string: str, fname: path_type, options: Optional[dict] = None) -> Union[str, None]:
    """
    Wrapper for command-line grep.

    Can pass any grep options like: options = {'A': value}

    :param str string: Search string
    :param str fname: File name to search
    :param Optional[dict] options: Grep options.
    :return output: String if matched, None if failed.
    """
    opts = ""
    if options is not None:
        for key, value in options.items():
            opts += "-" + key + " " + str(value) + " "

    grep_str = "grep " + opts + " '" + string + "' " + str(fname)

    try:
        output = subprocess.check_output(grep_str, shell=True).decode("utf-8")
    except subprocess.CalledProcessError as grepexc:
        warnings.warn(f"subprocess error: {grepexc.returncode}, grep found: {grepexc.output}")
        output = None

    return output
