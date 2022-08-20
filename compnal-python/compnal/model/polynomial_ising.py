from typing import Union
from base_compnal import base_model
from compnal.lattice.one_dimension import Chain
from compnal.lattice.two_dimension import Square, Triangle, Honeycomb
from compnal.lattice.three_dimension import Cubic
from compnal.lattice.higher_dimension import InfiniteRange, AnyLattice

LatticeType = Union[Chain, Square, Triangle, Honeycomb, Cubic, InfiniteRange, AnyLattice] 

class PolynomialIsing():
    
    def __init__(self, lattice: LatticeType, interaction: dict[int, float]) -> None:
        self.__base_model = base_model.make_polynomial_ising(lattice=lattice, interaction=interaction)

    


