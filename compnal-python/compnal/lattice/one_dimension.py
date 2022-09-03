from base_compnal import base_lattice
from compnal.lattice.boundary_condition import (
    BoundaryCondition,
    cast_base_boundary_condition,
    cast_boundary_condition,
)


class Chain(base_lattice.Chain):
    """One-dimensional Chain lattice.

    Attributes:
        system_size (int): The system size.
        boundary_condition (BoundaryCondition): The boundary condition.
    """

    def __init__(
        self,
        system_size: int,
        boundary_condition: BoundaryCondition = BoundaryCondition.OBC,
    ) -> None:
        """
        Args:
            system_size (int): The system size.
            boundary_condition (BoundaryCondition, optional): The boundary condition. Defaults to BoundaryCondition.OBC.
        """
        super().__init__(
            system_size=system_size,
            boundary_condition=cast_boundary_condition(boundary_condition),
        )

    def get_system_size(self) -> int:
        """Get the system size.

        Returns:
            int: The system size.
        """
        return super().get_system_size()

    def set_system_size(self, system_size: int) -> None:
        """Set the system size.

        Args:
            system_size (int): The system size.
        """
        super().set_system_size(system_size)

    def get_boundary_condition(self) -> BoundaryCondition:
        """Get the boundary condition.

        Returns:
            BoundaryCondition: The boundary condition.
        """
        return cast_base_boundary_condition(super().get_boundary_condition())

    def set_boundary_condition(self, boundary_condition: BoundaryCondition) -> None:
        """Set the boundary condition.

        Args:
            boundary_condition (BoundaryCondition): The boundary condition.
        """
        super().set_boundary_condition(cast_boundary_condition(boundary_condition))

    @property
    def system_size(self) -> int:
        return self.get_system_size()

    @system_size.setter
    def system_size(self, system_size: int) -> None:
        self.set_system_size(system_size)

    @property
    def boundary_condition(self) -> BoundaryCondition:
        return self.get_boundary_condition()

    @boundary_condition.setter
    def boundary_condition(self, bc: BoundaryCondition) -> None:
        self.set_boundary_condition(bc)
