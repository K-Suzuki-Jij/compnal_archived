from enum import Enum
from base_compnal import base_solver

class Updater(Enum):
    METROPOLIS = 0
    HEAT_BATH = 1

def cast_base_updater(updater: base_solver.CMCUpdater):
    if updater == base_solver.CMCUpdater.METROPOLIS:
        return Updater.METROPOLIS
    elif updater == base_solver.CMCUpdater.HEAT_BATH:
        return Updater.HEAT_BATH
    else:
        raise RuntimeError("Unknonw updater")

def cast_updater(updater: Updater):
    if updater == Updater.METROPOLIS:
        return base_solver.CMCUpdater.METROPOLIS
    elif updater == Updater.HEAT_BATH:
        return base_solver.CMCUpdater.HEAT_BATH
    else:
        raise RuntimeError("Unknonw updater")