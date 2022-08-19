from base_compnal import base_lattice
from compnal.lattice.boundary_condition import (
    BoundaryCondition, 
    cast_boundary_condition,
    cast_base_boundary_condition
)

class Chain(base_lattice.Chain):

    def __init__(
        self, 
        system_size: int, 
        boundary_condition: BoundaryCondition = BoundaryCondition.OBC
    ) -> None:

        super().__init__(
            system_size=system_size, 
            boundary_condition=cast_boundary_condition(boundary_condition)
        )

    @property
    def system_size(self):
        return super().get_system_size()

    @system_size.setter
    def system_size(self, system_size):
        super().set_system_size(system_size)

    @property
    def boundary_condition(self):
        return cast_base_boundary_condition(
            super().get_boundary_condition()
        )
        

    @boundary_condition.setter
    def boundary_condition(self, bc: BoundaryCondition):
        super().set_boundary_condition(
            cast_boundary_condition(bc)
        )