from typing import Union

from base_compnal import base_model
from compnal.lattice.boundary_condition import (
    BoundaryCondition,
    cast_base_boundary_condition,
)
from compnal.lattice.one_dimension import Chain
from compnal.lattice.two_dimension import Honeycomb, Square, Triangle
from compnal.lattice.three_dimension import Cubic
from compnal.lattice.higher_dimension import AnyLattice, InfiniteRange


class Ising:
    """Ising class.

    Attributes:
        system_size (int): The system size
        boundary_condition (BoundaryCondition): The boundary condition.
    """

    def __init__(
        self,
        lattice: Union[Chain, Square, Triangle, Honeycomb, Cubic, InfiniteRange],
        linear: float,
        quadratic: float,
    ) -> None:
        """Constructor.

        Args:
            lattice (Union[Chain, Square, Triangle, Honeycomb, Cubic, InfiniteRange]): The lattice.
            linear (float): The linear interaction.
            quadratic (float): The quadratic interaction.
        """

        self.__base_model = base_model.make_ising(
            lattice=lattice, interaction_deg_1=linear, interaction_deg_2=quadratic
        )

    def set_constant(self, constant: float) -> None:
        """Set the constant term.

        Args:
            constant (float): The constant term.
        """
        self.__base_model.set_constant(constant)

    def get_interaction(self) -> tuple[float, float, float]:
        """Get interaction.

        Returns:
            list:
                The interaction as tuple of the first value being constant term,
                the second value being linear term, and the third value being quadratic term.
        """
        return self.__base_model.get_interaction()

    def get_system_size(self) -> int:
        """Get the system size.

        Returns:
            int: The system size.
        """
        return self.__base_model.get_system_size()

    def get_boundary_condition(self) -> BoundaryCondition:
        """Get the boundary condition.

        Returns:
            BoundaryCondition: The boundary condition.
        """
        return cast_base_boundary_condition(self.__base_model.get_boundary_condition())

    def get_degree(self) -> int:
        """Get the degree of the interactions.

        Returns:
            int: The degree of the interactions.
        """
        return self.__base_model.get_degree()

    def calculate_energy(self, sample: list[int]) -> float:
        """Calculate energy.

        Args:
            sample (list[int]): Spin configuration.

        Returns:
            float: The energy.
        """
        return self.__base_model.calculate_energy(sample)

    @property
    def system_size(self) -> int:
        return self.get_system_size()

    @property
    def boundary_condition(self) -> BoundaryCondition:
        return self.get_boundary_condition()

    @property
    def _base_model(self) -> None:
        return self.__base_model


class IsingAnyLattice:
    """Ising class with the any lattice.

    Attributes:
        system_size (int): The system size
        boundary_condition (BoundaryCondition): The boundary condition.
    """

    def __init__(
        self,
        lattice: AnyLattice,
        linear: dict[Union[int, str, tuple[Union[int, str]]], float],
        quadratic: dict[
            tuple[
                Union[int, str, tuple[Union[int, str]]],
                Union[int, str, tuple[Union[int, str]]],
            ],
            float,
        ],
    ) -> None:
        """Constructor

        Args:
            lattice (AnyLattice): The lattice.
            linear (dict[Union[int, str, tuple[Union[int, str]]], float]): The linear interaction.
            quadratic (dict[tuple[Union[int, str, tuple[Union[int, str]]], Union[int, str, tuple[Union[int, str]]]], float]): The quadratic interaction.
        """
        self.__base_model = base_model.make_ising(
            lattice=lattice, linear=linear, quadratic=quadratic
        )

    def set_constant(self, constant: float) -> None:
        """Set the constant term.

        Args:
            constant (float): The constant term.
        """
        self.__base_model.set_constant(constant)

    def get_interaction(
        self,
    ) -> dict[tuple[Union[int, str, tuple[Union[int, str]]]], float]:
        """Get interaction.

        Returns:
            dict[tuple[Union[int, str, tuple[Union[int, str]]]], float]: The interaction.
        """
        interaction_map = {}
        interaction_map[()] = self.__base_model.get_constant()

        keys, values = self.__base_model.generate_linear_interaction_as_pair()
        for key, value in zip(keys, values):
            if isinstance(key, list):
                interaction_map[(tuple(key),)] = value
            else:
                interaction_map[(key,)] = value

        keys, values = self.__base_model.generate_quadratic_interaction_as_pair()
        for key, value in zip(keys, values):
            key_1 = tuple(key[0]) if isinstance(key[0], list) else key[0]
            key_2 = tuple(key[1]) if isinstance(key[1], list) else key[1]
            interaction_map[(key_1, key_2)] = value

        return interaction_map

    def get_system_size(self) -> int:
        """Get the system size.

        Returns:
            int: The system size.
        """
        return self.__base_model.get_system_size()

    def get_boundary_condition(self) -> BoundaryCondition:
        """Get the boundary condition.

        Returns:
            BoundaryCondition: The boundary condition.
        """
        return cast_base_boundary_condition(self.__base_model.get_boundary_condition())

    def get_degree(self) -> int:
        """Get the degree of the interactions.

        Returns:
            int: The degree of the interactions.
        """
        return self.__base_model.get_degree()

    def calculate_energy(self, sample: list[int]) -> float:
        """Calculate energy.

        Args:
            sample (list[int]): Spin configuration.

        Returns:
            float: The energy.
        """
        return self.__base_model.calculate_energy(sample)

    def generate_index_list(self) -> list[Union[int, str, list[Union[int, str]]]]:
        """Generate index list.

        Returns:
            list[Union[int, str, list[Union[int, str]]]]: The index list.
        """
        return self.__base_model.generate_index_list()

    @property
    def system_size(self) -> int:
        return self.get_system_size()

    @property
    def boundary_condition(self) -> BoundaryCondition:
        return self.get_boundary_condition()

    @property
    def _base_model(self):
        return self.__base_model


def make_ising(
    lattice: Union[
        Chain, Square, Triangle, Honeycomb, Cubic, InfiniteRange, AnyLattice
    ],
    linear: Union[float, dict[Union[int, str, list[Union[int, str]]], float]],
    quadratic: Union[
        float,
        dict[
            tuple[
                Union[int, str, list[Union[int, str]]],
                Union[int, str, list[Union[int, str]]],
            ],
            float,
        ],
    ],
) -> Union[Ising, IsingAnyLattice]:
    """Make PolynomialIsing class from lattice (and interaction).

    Args:
        lattice (Union[Chain, Square, Triangle, Honeycomb, Cubic, InfiniteRange, AnyLattice]):
            The lattice.
        linear (Union[float, dict[Union[int, str, list[Union[int, str]]], float]]):
            The linear interaction.
        quadratic (Union[float, dict[tuple[Union[int, str, list[Union[int, str]]], Union[int, str, list[Union[int, str]]]], float]]):
            The quadratic interaction.

    Raises:
        TypeError: The error raises when the unsupported lattices input.

    Returns:
        Union[Ising, IsingAnyLattice]: Ising class.
    """

    if isinstance(lattice, (Chain, Square, Triangle, Honeycomb, Cubic, InfiniteRange)):
        return Ising(lattice=lattice, linear=linear, quadratic=quadratic)
    elif isinstance(lattice, AnyLattice):
        return IsingAnyLattice(lattice=lattice, linear=linear, quadratic=quadratic)
    else:
        raise TypeError("Unknown LatticeType.")
