from typing import Optional, Union
from enum import Enum
from base_compnal import base_solver
from compnal.model.polynomial_ising import PolynomialIsing, PolynomialIsingAnyLattice
from compnal.solver.updater import Updater, cast_updater, cast_base_updater

ModelType = Union[PolynomialIsing, PolynomialIsingAnyLattice]

class ClassicalMonteCarlo:

    def __init__(
            self,
            model: ModelType,
            updater: Updater = Updater.METROPOLIS
        ) -> None:
        
        self.__base_solver = base_solver.make_classical_monte_carlo(
            model=model._base_model, 
            cmc_updater=cast_updater(updater)
        )
        
    def set_num_sweeps(self, num_sweeps: int) -> None:
        self.__base_solver.set_num_sweeps(num_sweeps=num_sweeps)

    def set_num_samples(self, num_samples: int) -> None:
        self.__base_solver.set_num_samples(num_samples=num_samples)

    def set_temperature(self, temperature: int) -> None:
        self.__base_solver.set_temperature(temperature=temperature)

    def set_inverse_temperature(self, beta: float) -> None:
        self.__base_solver.set_inverse_temperature(inverse_temperature=beta)
    
    def get_num_sweeps(self) -> int:
        return self.__base_solver.get_num_sweeps()

    def get_num_samples(self) -> int:
        return self.__base_solver.get_num_samples()

    def get_temperature(self) -> int:
        return self.__base_solver.get_temperature()

    def get_inverse_temperature(self) -> float:
        return self.__base_solver.get_inverse_temperature()

    def get_seed(self) -> int:
        return self.__base_solver.get_seed()

    def run(self, seed: Optional[int] = None) -> int:
        if seed is None:
            self.__base_solver.run()
        else:
            self.__base_solver.run(seed=seed)

    def calculate_sample_average(self) -> float:
        return self.__base_solver.calculate_sample_average()

    def calculate_sample_moment(self, degree: int) -> float:
        return self.__base_solver.calculate_sample_moment(degree=degree)

    @property
    def updater(self) -> Updater:
        return cast_base_updater(self.__base_solver.cmc_updater)

    @updater.setter
    def updater(self, updater: Updater) -> None:
        self.__base_solver.cmc_updater = cast_updater(updater)

    @property
    def num_sweeps(self) -> int:
        return self.get_num_sweeps()

    @num_sweeps.setter
    def num_sweeps(self, num_sweeps: int) -> None:
        self.set_num_sweeps(num_sweeps=num_sweeps)

    @property
    def num_samples(self) -> int:
        return self.get_num_samples()

    @num_samples.setter
    def num_samples(self, num_samples: int) -> None:
        self.set_num_samples(num_samples=num_samples)

    @property
    def beta(self) -> float:
        return self.get_inverse_temperature()

    @beta.setter
    def beta(self, beta: float) -> None:
        self.set_inverse_temperature(beta=beta)
    