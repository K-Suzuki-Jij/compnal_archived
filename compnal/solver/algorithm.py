from enum import Enum

from compnal.base_compnal import base_solver


class Algorithm(Enum):
    """Update algorithm."""

    METROPOLIS = 0
    HEAT_BATH = 1


def cast_base_algorithm(algorithm: base_solver.Algorithm) -> Algorithm:
    """Cast the algorithm from base_solver.Algorithm to solver.Algorithm.

    Args:
        algorithm (base_solver.Algorithm): The algorithm as base_solver.Algorithm.

    Raises:
        RuntimeError: When the unknown algorithm is input.

    Returns:
        Algorithm: The algorithm as solver.Algorithm.
    """
    if algorithm == base_solver.Algorithm.METROPOLIS:
        return Algorithm.METROPOLIS
    elif algorithm == base_solver.Algorithm.HEAT_BATH:
        return Algorithm.HEAT_BATH
    else:
        raise RuntimeError("Unknown algorithm")


def cast_algorithm(algorithm: Algorithm):
    """Cast the algorithm from solver.Algorithm to base_solver.Algorithm.

    Args:
        algorithm (Algorithm): The algorithm as solver.Algorithm.

    Raises:
        RuntimeError: When the unknown algorithm is input.

    Returns:
        base_solver.Algorithm: The algorithm as base_solver.Algorithm.
    """
    if algorithm == Algorithm.METROPOLIS:
        return base_solver.Algorithm.METROPOLIS
    elif algorithm == Algorithm.HEAT_BATH:
        return base_solver.Algorithm.HEAT_BATH
    else:
        raise RuntimeError("Unknown algorithm")
