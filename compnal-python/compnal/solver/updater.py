from enum import Enum

from base_compnal import base_solver


class Updater(Enum):
    """Update algorithm."""

    METROPOLIS = 0
    HEAT_BATH = 1


def cast_base_updater(updater: base_solver.CMCUpdater) -> Updater:
    """Cast the updater from base_solver.CMCUpdater to solver.Updater.

    Args:
        updater (base_solver.CMCUpdater): The updater as base_solver.CMCUpdater.

    Raises:
        RuntimeError: When the unknown updater is input.

    Returns:
        Updater: The updater as solver.Updater.
    """
    if updater == base_solver.CMCUpdater.METROPOLIS:
        return Updater.METROPOLIS
    elif updater == base_solver.CMCUpdater.HEAT_BATH:
        return Updater.HEAT_BATH
    else:
        raise RuntimeError("Unknown updater")


def cast_updater(updater: Updater):
    """Cast the updater from solver.Updater to base_solver.CMCUpdater.

    Args:
        updater (Updater): The updater as solver.Updater.

    Raises:
        RuntimeError: When the unknown updater is input.

    Returns:
        base_solver.CMCUpdater: The updater as base_solver.CMCUpdater.
    """
    if updater == Updater.METROPOLIS:
        return base_solver.CMCUpdater.METROPOLIS
    elif updater == Updater.HEAT_BATH:
        return base_solver.CMCUpdater.HEAT_BATH
    else:
        raise RuntimeError("Unknown updater")
