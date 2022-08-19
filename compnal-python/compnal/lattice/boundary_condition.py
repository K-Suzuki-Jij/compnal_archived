from enum import Enum
from base_compnal import base_lattice

class BoundaryCondition(Enum):
    NONE = 0
    OBC = 1
    PBC = 2

def cast_base_boundary_condition(boundary_condition: base_lattice.BoundaryCondition):
    if boundary_condition == base_lattice.BoundaryCondition.NONE:
        return BoundaryCondition.NONE
    elif boundary_condition == base_lattice.BoundaryCondition.OBC:
        return BoundaryCondition.OBC
    elif boundary_condition == base_lattice.BoundaryCondition.PBC:
        return BoundaryCondition.PBC
    else:
        raise RuntimeError("Unknonw BoundaryCondition")

def cast_boundary_condition(boundary_condition: BoundaryCondition):
    if boundary_condition == BoundaryCondition.NONE:
        return base_lattice.BoundaryCondition.NONE
    elif boundary_condition == BoundaryCondition.OBC:
        return base_lattice.BoundaryCondition.OBC
    elif boundary_condition == BoundaryCondition.PBC:
        return base_lattice.BoundaryCondition.PBC
    else:
        raise RuntimeError("Unknonw BoundaryCondition")



