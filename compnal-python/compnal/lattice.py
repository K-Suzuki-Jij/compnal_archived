from base_compnal import lattice
from base_compnal.lattice import BoundaryCondition as BC

class Chain(lattice.Chain):

    def __init__(self, system_size: int, bc: BC = BC.OBC) -> None:
        super().__init__(system_size=system_size, boundary_condition=bc)

    @property
    def system_size(self):
        return super().get_system_size()

    @property
    def boundary_condition(self):
        return super().get_boundary_condition()
    

class Square(lattice.Square):
    pass




