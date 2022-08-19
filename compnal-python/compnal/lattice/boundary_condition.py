from enum import Enum
from base_compnal import base_lattice

class BoundaryCondition(Enum):
    """Boundary condition.

    Args:
        NONE: None type.
        OBC: Open boundary condition.
        PBC: Periodic boundary condition.
    """
    NONE = 0
    OBC = 1
    PBC = 2

def cast_base_boundary_condition(boundary_condition: base_lattice.BoundaryCondition):
    """_summary_

    Args:
        boundary_condition (base_lattice.BoundaryCondition): _description_

    Raises:
        RuntimeError: _description_

    Returns:
        _type_: _description_
    """
    if boundary_condition == base_lattice.BoundaryCondition.NONE:
        return BoundaryCondition.NONE
    elif boundary_condition == base_lattice.BoundaryCondition.OBC:
        return BoundaryCondition.OBC
    elif boundary_condition == base_lattice.BoundaryCondition.PBC:
        return BoundaryCondition.PBC
    else:
        raise RuntimeError("Unknonw BoundaryCondition")

def cast_boundary_condition(boundary_condition: BoundaryCondition):
    """_summary_

    Args:
        boundary_condition (BoundaryCondition): _description_

    Raises:
        RuntimeError: _description_

    Returns:
        _type_: _description_
    """
    if boundary_condition == BoundaryCondition.NONE:
        return base_lattice.BoundaryCondition.NONE
    elif boundary_condition == BoundaryCondition.OBC:
        return base_lattice.BoundaryCondition.OBC
    elif boundary_condition == BoundaryCondition.PBC:
        return base_lattice.BoundaryCondition.PBC
    else:
        raise RuntimeError("Unknonw BoundaryCondition")



