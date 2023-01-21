"""Utils for serialization for jobflow.
"""
import os


def special_serialization_attrs(instance) -> dict:
    """Gives special keys + values for serialization.
     Currently only supports jobflow via the env_var USE_JOBFLOW.

    :param instance: object you want to serialise
    :returns: dictionary with the special keys + values
    """
    if os.getenv('USE_JOBFLOW') is not None:
        return {"@module": instance.__class__.__module__, "@class": instance.__class__.__name__, }
    return {}
