"""Base class for inheriting for serialisable objects."""

import os

try:
    from monty.json import MSONable

    supercls = MSONable
except ModuleNotFoundError:
    supercls = object

# not sure at the moment why anyone would like to do this
# if this changes after importing excitingtools, the module must be reloaded
if os.getenv("USE_MONTY") == "false":
    supercls = object


class ECTObject(supercls):
    """Base excitingtools object."""
