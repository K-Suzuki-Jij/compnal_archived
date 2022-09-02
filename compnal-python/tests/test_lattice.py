from compnal import lattice

def test_lattice_chain():
    chain = lattice.Chain(system_size=10, boundary_condition=lattice.BoundaryCondition.OBC)
    chain.system_size = 7
    chain.boundary_condition = lattice.BoundaryCondition.PBC
    assert chain.get_system_size() == 7
    assert chain.system_size == 7
    assert chain.get_boundary_condition() == lattice.BoundaryCondition.PBC
    assert chain.boundary_condition == lattice.BoundaryCondition.PBC

def test_lattice_square():
    square = lattice.Square(x_size=3, y_size=4, boundary_condition=lattice.BoundaryCondition.PBC)
    square.x_size = 2
    square.y_size = 3
    square.boundary_condition = lattice.BoundaryCondition.OBC
    assert square.get_system_size() == 6
    assert square.system_size == 6
    assert square.get_boundary_condition() == lattice.BoundaryCondition.OBC
    assert square.boundary_condition == lattice.BoundaryCondition.OBC