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

        return PolynomialIsing(lattice=lattice, interaction=interaction)
                

    elif isinstance(lattice, AnyLattice):
        class PolynomialIsing(base_poly_ising_class):

            def __init__(
                self, 
                lattice: LatticeType, 
                interaction: Dict[List[Union[int, str, List[Union[int, str]]]], float] = {}
            ) -> None:

                super().__init__(lattice=lattice, interaction=interaction)

        return PolynomialIsing(lattice=lattice, interaction=interaction)
    else:
        raise TypeError("Unknown LatticeType.")



    


