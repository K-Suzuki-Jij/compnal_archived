from typing import Union, Dict, List
from base_compnal import base_model
from compnal.lattice.one_dimension import Chain
from compnal.lattice.two_dimension import Square, Triangle, Honeycomb
from compnal.lattice.three_dimension import Cubic
from compnal.lattice.higher_dimension import InfiniteRange, AnyLattice

BasicLatticeType = Union[Chain, Square, Triangle, Honeycomb, Cubic, InfiniteRange]
LatticeType = Union[BasicLatticeType, AnyLattice]
InteractionType = Union[Dict[int, float], Dict[List[Union[int, str, List[Union[int, str]]]], float]]


def make_polynomial_ising(lattice: LatticeType, interaction: InteractionType = {}):
    
    base_poly_ising_class = base_model.make_polynomial_ising(lattice=lattice).__class__

    if isinstance(lattice, (Chain, Square, Triangle, Honeycomb, Cubic, InfiniteRange)):

        class PolynomialIsing(base_poly_ising_class):

            def __init__(
                self, 
                lattice: LatticeType, 
                interaction: Dict[int, float] = {}
            ) -> None:

                super().__init__(lattice=lattice, interaction=interaction)

            def set_interaction(self, degree: int, value: float) -> None:
                super().set_interaction(degree, value)

            def add_interaction(self, degree: int, value: float) -> None:
                super().add_interaction(degree, value)

            def get_interaction(self) -> list:
                return super().get_interaction()

            def get_system_size(self) -> int:
                return super().get_system_size()

            def calculate_energy(self, sample: list[int]) -> float:
                return super().calculate_energy(sample)

            @property
            def system_size(self) -> int:
                return self.get_system_size()
            
        return PolynomialIsing(lattice=lattice, interaction=interaction)
                

    elif isinstance(lattice, AnyLattice):
        class PolynomialIsing(base_poly_ising_class):

            def __init__(
                self, 
                lattice: AnyLattice, 
                interaction: Dict[List[Union[int, str, List[Union[int, str]]]], float] = {}
            ) -> None:

                super().__init__(lattice=lattice, interaction=interaction)

            def set_interaction(self, index_list: List[Union[int, str, List[Union[int, str]]]], value: float) -> None:
                super().set_interaction(index_list, value)

            def add_interaction(self, index_list: List[Union[int, str, List[Union[int, str]]]], value: float) -> None:
                super().add_interaction(index_list, value)

            def get_interaction(self) -> list:
                keys, values = super().generate_interaction_as_pair()
                interaction_map = {}
                for key, value in zip(keys, values):
                    if isinstance(key, list):
                        interaction_map[tuple(key)] = value
                    else:
                        interaction_map[key] = value
                return interaction_map

            def get_system_size(self) -> int:
                return super().get_system_size()

            def calculate_energy(self, sample: list[int]) -> float:
                return super().calculate_energy(sample)

            def generate_index_list(self) -> List[Union[int, str, List[Union[int, str]]]]:
                return super().generate_index_list()

            @property
            def system_size(self) -> int:
                return self.get_system_size()

        return PolynomialIsing(lattice=lattice, interaction=interaction)
    else:
        raise TypeError("Unknown LatticeType.")



    


