from typing import Union

from base_compnal import base_model
from compnal.lattice.boundary_condition import (
    BoundaryCondition,
    cast_base_boundary_condition,
)
from compnal.lattice.higher_dimension import AnyLattice, InfiniteRange
from compnal.lattice.one_dimension import Chain
from compnal.lattice.three_dimension import Cubic
from compnal.lattice.two_dimension import Honeycomb, Square, Triangle

LatticeType = Union[
    Chain, Square, Triangle, Honeycomb, Cubic, InfiniteRange, AnyLattice
]
InteractionType = Union[
    dict[int, float], dict[tuple[Union[int, str, tuple[Union[int, str]]]], float]
]


class PolynomialIsing:
    """PolynomialIsing class.

    Attributes:
        system_size (int): The system size
        boundary_condition (BoundaryCondition): The boundary condition.
    """

    def __init__(
        self,
        lattice: Union[Chain, Square, Triangle, Honeycomb, Cubic, InfiniteRange],
        interaction: dict[int, float],
    ) -> None:
        """
        Args:
            lattice (Union[Chain, Square, Triangle, Honeycomb, Cubic, InfiniteRange]): The lattice.
            interaction (dict[int, float]): The interaction.
        """
        self.__base_model = base_model.make_polynomial_ising(
            lattice=lattice, interaction=interaction
        )

    def get_interaction(self) -> list:
        """Get interaction.

        Returns:
            list: The interaction.
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
    def _base_model(self):
        return self.__base_model


class PolynomialIsingAnyLattice:
    """PolynomialIsing class with the any lattice.
    User can add any interactions.
    """

    def __init__(
        self,
        lattice: AnyLattice,
        interaction: dict[tuple[Union[int, str, tuple[Union[int, str]]]], float],
    ) -> None:
        """Constructor of PolynomialIsingAnyLattice class.

        Args:
            lattice (AnyLattice): The lattice.
            interaction (dict[tuple[Union[int, str, tuple[Union[int, str]]]], float]): The interactions.
        """

        self.__base_model = base_model.make_polynomial_ising(
            lattice=lattice, interaction=interaction
        )

    def get_interaction(
        self,
    ) -> dict[tuple[Union[int, str, tuple[Union[int, str]]]], float]:
        """Get interaction.

        Returns:
            dict[tuple[Union[int, str, tuple[Union[int, str]]]], float]: The interaction.
        """
        keys, values = self.__base_model.generate_interaction_as_pair()
        interaction_map = {}
        for key, value in zip(keys, values):
            if isinstance(key, list):
                interaction_map[tuple(key)] = value
            else:
                interaction_map[key] = value
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

    def get_index_list(self) -> list[Union[int, str, list[Union[int, str]]]]:
        """Generate index list.

        Returns:
            list[Union[int, str, list[Union[int, str]]]]: The index list.
        """
        return self.__base_model.get_index_list()

    @property
    def system_size(self) -> int:
        return self.__base_model.get_system_size()

    @property
    def boundary_condition(self) -> BoundaryCondition:
        return self.__base_model.get_boundary_condition()

    @property
    def _base_model(self):
        return self.__base_model


def make_polynomial_ising(
    lattice: LatticeType, interaction: InteractionType
) -> Union[PolynomialIsing, PolynomialIsingAnyLattice]:
    """Make PolynomialIsing class from lattice (and interaction).

    Args:
        lattice (LatticeType): The lattice.
        interaction (InteractionType): The interaction.

    Raises:
        TypeError: The error raises when the unsupported lattices input.

    Returns:
        Union[PolynomialIsing, PolynomialIsingAnyLattice]: PolynomialIsing class.
    """

    if isinstance(lattice, (Chain, Square, Triangle, Honeycomb, Cubic, InfiniteRange)):
        return PolynomialIsing(lattice=lattice, interaction=interaction)
    elif isinstance(lattice, AnyLattice):
        return PolynomialIsingAnyLattice(lattice=lattice, interaction=interaction)
    else:
        raise TypeError("Unknown LatticeType.")
