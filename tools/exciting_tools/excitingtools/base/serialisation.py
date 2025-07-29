"""Base class for inheriting for serialisable objects and models.

Note that I left out pydantic as an optional dependency in the pyproject.toml file. It's not worth adding it as a
separate field as you can just install it manually. Anyway, the intended use case is together with excitingworkflow
and jobflow(-remote), and there you have it installed for sure.
"""

import os
from importlib.util import find_spec

pydantic_avail = find_spec("pydantic") is not None
monty_avail = find_spec("monty") is not None

if monty_avail:
    from monty.json import MSONable

    ECTModelBase, ECTObject = MSONable, MSONable
else:
    ECTModelBase, ECTObject = object, object  # Just a dummy base class

if pydantic_avail:
    from pydantic import BaseModel

    ECTModelBase = BaseModel

    def model_decorator(cls):
        return cls  # No-op: Pydantic doesn't need a decorator
else:
    from dataclasses import dataclass

    def model_decorator(cls):
        return dataclass(cls)  # Apply @dataclass dynamically


# not sure at the moment why anyone would like to do this
# if this changes after importing excitingtools, the module must be reloaded
if os.getenv("USE_MONTY") == "false":
    ECTObject = object

__all__ = ["ECTModelBase", "ECTObject", "model_decorator"]
