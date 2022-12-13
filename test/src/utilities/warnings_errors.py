""" Utilities for raising warnings and errors.
"""
from typing import Type
import warnings


def warnings_as_errors(msg: str, exception: Type[Exception], is_error=False):
    """Raise a warning or error.

    Call as warnings_as_errors(msg, ValueError, settings.SAFE_MODE)

    Note, could have done something like:
    ```
       def warnings_as_errors(msg: str, error: False):
       error_value = {True: "error", False: "default"}
       with warnings.catch_warnings():
           warnings.simplefilter(error_value[error])
           warnings.warn(message=msg)
    ```
    but this is more code.

    :param msg: User message
    :param exception: Exception type (not instance) i.e. ValueError is the type
    ValueError() is an instance of that type
    :param is_error: Treat the warning as an error
    """
    if is_error:
        raise exception(msg)
    else:
        warnings.warn(message=msg)
