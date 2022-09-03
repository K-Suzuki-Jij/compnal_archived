from compnal import lattice, model


def test_model_make_polynomial_ising():
    ising = model.make_polynomial_ising(
        lattice=lattice.Chain(system_size=10), interaction={1: 2.0, 3: 7.1}
    )
    assert isinstance(ising, model.PolynomialIsing)

    ising = model.make_polynomial_ising(
        lattice=lattice.Square(x_size=3, y_size=4), interaction={1: 2.0, 3: 7.1}
    )
    assert isinstance(ising, model.PolynomialIsing)

    ising = model.make_polynomial_ising(
        lattice=lattice.Triangle(x_size=3, y_size=4), interaction={1: 2.0, 3: 7.1}
    )
    assert isinstance(ising, model.PolynomialIsing)

    ising = model.make_polynomial_ising(
        lattice=lattice.Honeycomb(x_size=3, y_size=4), interaction={1: 2.0, 3: 7.1}
    )
    assert isinstance(ising, model.PolynomialIsing)

    ising = model.make_polynomial_ising(
        lattice=lattice.Cubic(x_size=3, y_size=4, z_size=2),
        interaction={1: 2.0, 3: 7.1},
    )
    assert isinstance(ising, model.PolynomialIsing)

    ising = model.make_polynomial_ising(
        lattice=lattice.InfiniteRange(system_size=4), interaction={1: 2.0, 3: 7.1}
    )
    assert isinstance(ising, model.PolynomialIsing)

    ising = model.make_polynomial_ising(
        lattice=lattice.AnyLattice(),
        interaction={(1, "a"): -2.0, (1, (2, "a"), "b"): -3.1},
    )
    assert isinstance(ising, model.PolynomialIsingAnyLattice)
