from base_compnal import base_lattice
from compnal.lattice.boundary_condition import (
    BoundaryCondition,
    cast_base_boundary_condition,
    cast_boundary_condition,
)


class Cubic(base_lattice.Cubic):
    """Three-dimensional cubic lattice.

    Attributes:
        x_size (int): The size of the x-direction.
        y_size (int): The size of the y-direction.
        z_size (int): The size of the z-direction.
        system_size (int): The system size (x_size*y_size*z_size).
        boundary_condition (BoundaryCondition): The boundary condition.
    """

    def __init__(
        self,
        x_size: int,
        y_size: int,
        z_size: int,
        boundary_condition: BoundaryCondition = BoundaryCondition.OBC,
    ) -> None:
        """
        Args:
            x_size (int): The size of the x-direction.
            y_size (int): The size of the y-direction.
            z_size (int): The size of the z-direction.
            boundary_condition (BoundaryCondition, optional): The boundary condition. Defaults to BoundaryCondition.OBC.
        """
        super().__init__(
            x_size=x_size,
            y_size=y_size,
            z_size=z_size,
            boundary_condition=cast_boundary_condition(boundary_condition),
        )

    def get_x_size(self) -> int:
        """Get the size of the x-direction.

        Returns:
            int: The size of the x-direction.
        """
        return super().get_x_size()

    def get_y_size(self) -> int:
        """Get the size of the y-direction.

        Returns:
            int: The size of the y-direction.
        """
        return super().get_y_size()

    def get_z_size(self) -> int:
        """Get the size of the z-direction.

        Returns:
            int: The size of the z-direction.
        """
        return super().get_z_size()

    def get_system_size(self) -> int:
        """Get the system size.

        Returns:
            int: The system size.
        """
        return super().get_system_size()

    def get_boundary_condition(self) -> BoundaryCondition:
        """Get the boundary condition.

        Returns:
            BoundaryCondition: The boundary condition.
        """
        return cast_base_boundary_condition(super().get_boundary_condition())

    def set_x_size(self, x_size: int) -> None:
        """Set the size of the x-direction.

        Args:
            x_size (int): the size of the x-direction.

        Returns:
            None
        """
        super().set_x_size(x_size)

    def set_y_size(self, y_size: int) -> None:
        """Set the size of the y-direction.

        Args:
            y_size (int): the size of the y-direction.

        Returns:
            None
        """
        super().set_y_size(y_size)

    def set_z_size(self, z_size: int) -> None:
        """Set the size of the z-direction.

        Args:
            z_size (int): the size of the z-direction.

        Returns:
            None
        """
        super().set_z_size(z_size)

    def set_boundary_condition(self, boundary_condition: BoundaryCondition) -> None:
        """Set the boundary condition.

        Args:
            boundary_condition (BoundaryCondition): The boundary condition.

        Returns:
            None
        """
        super().set_boundary_condition(cast_boundary_condition(boundary_condition))

    @property
    def x_size(self) -> int:
        return self.get_x_size()

    @x_size.setter
    def x_size(self, x_size: int) -> None:
        self.set_x_size(x_size)

    @property
    def y_size(self) -> int:
        return self.get_y_size()

    @y_size.setter
    def y_size(self, y_size: int) -> None:
        self.set_y_size(y_size)

    @property
    def z_size(self) -> int:
        return self.get_z_size()

    @z_size.setter
    def z_size(self, z_size: int) -> None:
        self.set_z_size(z_size)

    @property
    def system_size(self) -> None:
        return self.get_system_size()

    @property
    def boundary_condition(self) -> BoundaryCondition:
        return self.get_boundary_condition()

    @boundary_condition.setter
    def boundary_condition(self, boundary_condition: BoundaryCondition) -> None:
        self.set_boundary_condition(boundary_condition)
