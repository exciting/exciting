"""Utils for serialization."""

import copy
import os


def special_serialization_attrs(instance) -> dict:
    """Gives special keys + values for serialization.

     Currently only supports monty/jobflow via the env_var USE_JOBFLOW.

    :param instance: object you want to serialize
    :returns: dictionary with the special keys + values
    """
    if os.getenv("USE_JOBFLOW") is not None:
        return {"@module": instance.__class__.__module__, "@class": instance.__class__.__name__}
    return {}


def deserialize_object(cls, d):
    """Recreates a class instance from a serialized dictionary.

    :param cls: the class
    :param d: the serialized dictionary
    :return: the class instance
    """
    my_dict = copy.deepcopy(d)
    # Remove key value pairs needed for workflow programs
    # call function on class to get only the keys (values not needed)
    serialise_keys = special_serialization_attrs(cls)
    for key in serialise_keys:
        my_dict.pop(key, None)
    return cls(**my_dict)
