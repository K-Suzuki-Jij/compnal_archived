from typing import Tuple, Union, Dict, List
from base_compnal import base_model
from compnal.lattice.one_dimension import Chain
from compnal.lattice.two_dimension import Square, Triangle, Honeycomb
from compnal.lattice.three_dimension import Cubic
from compnal.lattice.higher_dimension import InfiniteRange, AnyLattice
from compnal.lattice.boundary_condition import (
    BoundaryCondition, 
    cast_base_boundary_condition
)


BasicLatticeType = Union[Chain, Square, Triangle, Honeycomb, Cubic, InfiniteRange]
LatticeType = Union[BasicLatticeType, AnyLattice]
InteractionType = Union[Dict[int, float], Dict[List[Union[int, str, List[Union[int, str]]]], float]]


def make_polynomial_ising(lattice: LatticeType, interaction: InteractionType = {}):
    """Make PolynomialIsing class from lattice (and interaction).

    Args:
        lattice (LatticeType): The lattice.
        interaction (InteractionType, optional): The interaction. Defaults to {}.

    Raises:
        TypeError: The error raises when the unsupported lattices input.

    Returns:
        PolynomialIsing: PolynomialIsing class.
    """
    base_poly_ising_class = base_model.make_polynomial_ising(lattice=lattice).__class__

    if isinstance(lattice, (Chain, Square, Triangle, Honeycomb, Cubic, InfiniteRange)):

        class PolynomialIsing(base_poly_ising_class):
            """PolynomialIsing class.

            Attributes:
                system_size (int): The system size
                boundary_condition (BoundaryCondition): The boundary condition.
            """
            def __init__(
                self, 
                lattice: LatticeType, 
                interaction: Dict[int, float] = {}
            ) -> None:
                """
                Args:
                    lattice (LatticeType): The lattice.
                    interaction (Dict[int, float], optional): The interaction. Defaults to {}.
                """
                super().__init__(lattice=lattice, interaction=interaction)

            def set_interaction(self, degree: int, value: float) -> None:
                """Set interaction. Overwrite if already set.

                Args:
                    degree (int): The degree of the interaction.
                    value (float): The value of the interaction.
                """
                super().set_interaction(degree, value)

            def add_interaction(self, degree: int, value: float) -> None:
                """Add interaction. 

                Args:
                    degree (int): The degree of the interaction.
                    value (float): The value of the interaction.
                """
                super().add_interaction(degree, value)

            def get_interaction(self) -> list:
                """Get interaction.

                Returns:
                    list: The interaction.
                """
                return super().get_interaction()

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

            def calculate_energy(self, sample: list[int]) -> float:
                """Calculate energy.

                Args:
                    sample (list[int]): Spin configuration.

                Returns:
                    float: The energy.
                """
                return super().calculate_energy(sample)

            @property
            def system_size(self) -> int:
                return self.get_system_size()
            
            @property
            def boundary_condition(self) -> BoundaryCondition:
                return self.get_boundary_condition()
            
        return PolynomialIsing(lattice=lattice, interaction=interaction)
                

    elif isinstance(lattice, AnyLattice):
        class PolynomialIsing(base_poly_ising_class):
            """PolynomialIsing class. Users can set any interactions.

            Attributes:
                system_size (int): The system size
                boundary_condition (BoundaryCondition): The boundary condition.
            """
            def __init__(
                self, 
                lattice: AnyLattice, 
                interaction: Dict[List[Union[int, str, List[Union[int, str]]]], float] = {}
            ) -> None:
                """
                Args:
                    lattice (AnyLattice): The lattice.
                    interaction (Dict[List[Union[int, str, List[Union[int, str]]]], float]): The interaction. Defaults to {}.
                """
                super().__init__(lattice=lattice, interaction=interaction)

            def set_interaction(self, index_list: List[Union[int, str, List[Union[int, str]]]], value: float) -> None:
                """Set interaction. Overwrite if already set.

                Args:
                    index_list (List[Union[int, str, List[Union[int, str]]]]): The index list connected by the interaction.
                    value (float): The value of the interaction.
                """
                super().set_interaction(index_list, value)

            def add_interaction(self, index_list: List[Union[int, str, List[Union[int, str]]]], value: float) -> None:
                """Add interaction.

                Args:
                    index_list (List[Union[int, str, List[Union[int, str]]]]): The index list connected by the interaction.
                    value (float): The value of the interaction.
                """
                super().add_interaction(index_list, value)

            def get_interaction(self) -> Dict[List[Union[int, str, Tuple[Union[int, str]]]], float]:
                """Get interaction.

                Returns:
                    Dict[List[Union[int, str, Tuple[Union[int, str]]]], float]: The interaction.
                """
                keys, values = super().generate_interaction_as_pair()
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
                return super().get_system_size()

            def get_boundary_condition(self) -> BoundaryCondition:
                """Get the boundary condition.

                Returns:
                    BoundaryCondition: The boundary condition.
                """
                return cast_base_boundary_condition(super().get_boundary_condition())

            def calculate_energy(self, sample: list[int]) -> float:
                """Calculate energy.

                Args:
                    sample (list[int]): Spin configuration.

                Returns:
                    float: The energy.
                """
                return super().calculate_energy(sample)

            def generate_index_list(self) -> List[Union[int, str, List[Union[int, str]]]]:
                """Generate index list.

                Returns:
                    List[Union[int, str, List[Union[int, str]]]]: The index list.
                """
                return super().generate_index_list()

            @property
            def system_size(self) -> int:
                return self.get_system_size()

            @property
            def boundary_condition(self) -> BoundaryCondition:
                return self.get_boundary_condition()

        return PolynomialIsing(lattice=lattice, interaction=interaction)
    else:
        raise TypeError("Unknown LatticeType.")



    


