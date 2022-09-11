from typing import Optional, Union

from base_compnal import base_solver
from compnal.model.polynomial_ising import PolynomialIsing
from compnal.model.ising import Ising, IsingAnyLattice
from compnal.solver.updater import Updater, cast_base_updater, cast_updater

ModelType = Union[PolynomialIsing, Ising, IsingAnyLattice]


class ClassicalMonteCarlo:
    """ClassicalMonteCarlo class. This class can treat classical models.

    Attributes:
        updater (Updater): The algorithm of update method.
        num_sweeps (int): The number of sweeps.
        num_samples (int): The number of samples.
        num_threads (int): The number of threads.
        beta (float): The inverse temperature.
        model (Union[PolynomialIsing, Ising, IsingAnyLattice]): The model.
    """

    def __init__(self, model: ModelType, updater: Updater = Updater.METROPOLIS) -> None:
        """The constructor.

        Args:
            model (ModelType): The classical model.
            updater (Updater, optional): The algorithm of update method. Defaults to Updater.METROPOLIS.
        """

        self.__base_solver = base_solver.make_classical_monte_carlo(
            model=model._base_model, cmc_updater=cast_updater(updater)
        )

        self.__model = model

    def set_num_sweeps(self, num_sweeps: int) -> None:
        """Set the number of sweeps.

        Args:
            num_sweeps (int): The number of sweeps.
        """
        self.__base_solver.set_num_sweeps(num_sweeps=num_sweeps)

    def set_num_samples(self, num_samples: int) -> None:
        """Set the number of samples.

        Args:
            num_samples (int): The number of samples.
        """
        self.__base_solver.set_num_samples(num_samples=num_samples)

    def set_num_threads(self, num_threads: int) -> None:
        """Set the number of threads.

        Args:
            num_threads (int): The number of threads.
        """
        self.__base_solver.set_num_threads(num_threads=num_threads)

    def set_temperature(self, temperature: int) -> None:
        """Set the temperature at which the monte carlo algorithm executes.

        Args:
            temperature (int): The temperature.
        """
        self.__base_solver.set_temperature(temperature=temperature)

    def set_inverse_temperature(self, beta: float) -> None:
        """Set the inverse temperature at which the monte carlo algorithm executes.

        Args:
            beta (float): The inverse temperature.
        """
        self.__base_solver.set_inverse_temperature(inverse_temperature=beta)

    def get_num_sweeps(self) -> int:
        """Get the number of sweeps.

        Returns:
            int: The number of sweeps.
        """
        return self.__base_solver.get_num_sweeps()

    def get_num_samples(self) -> int:
        """Get the number of samples.

        Returns:
            int: The number of samples.
        """
        return self.__base_solver.get_num_samples()

    def get_num_threads(self) -> int:
        """Get the number of threads.

        Returns:
            int: The number of threads.
        """
        return self.__base_solver.get_num_threads()

    def get_temperature(self) -> int:
        """Get the temperature at which the monte carlo algorithm executes.

        Returns:
            int: The temperature.
        """
        return self.__base_solver.get_temperature()

    def get_inverse_temperature(self) -> float:
        """Get the inverse temperature at which the monte carlo algorithm executes.

        Returns:
            float: The inverse temperature.
        """
        return self.__base_solver.get_inverse_temperature()

    def get_seed(self) -> int:
        """Get the seed used by the monte carlo algorithm.

        Returns:
            int: The seed.
        """
        return self.__base_solver.get_seed()

    def get_samples(self) -> list[list[int]]:
        """Get the samples.

        Returns:
            list[list[int]]: The samples.
        """
        return self.__base_solver.get_samples()

    def run(self, seed: Optional[int] = None) -> None:
        """Run the monte carlo algorithm.

        Args:
            seed (Optional[int], optional): The seed used by the monte carlo algorithm. Defaults to None.
        """
        if seed is None:
            self.__base_solver.run()
        else:
            self.__base_solver.run(seed=seed)

    def calculate_sample_average(self) -> float:
        """Calculate average of all the samples.

        Returns:
            float: The average value.
        """
        return self.__base_solver.calculate_sample_average()

    def calculate_sample_distribution(self) -> list[float]:
        """Calculate average of samples at each sites.

        Returns:
            list[float]: Average of samples at each sites.
        """
        return self.__base_solver.calculate_sample_distribution()

    def calculate_sample_moment(self, degree: int) -> float:
        """Calculate the moment with the specified degree of each sample and take average for all the moments.

        Args:
            degree (int): The degree of the moment.

        Returns:
            float: The moment.
        """
        return self.__base_solver.calculate_sample_moment(degree=degree)

    def calculate_correlation(
        self,
        index_1: Union[int, list[Union[int, str, list[Union[int, str]]]]],
        index_2: Union[int, list[Union[int, str, list[Union[int, str]]]]],
    ) -> float:
        """Calculate correlation function.

        Args:
            index_1 (Union[int, list[Union[int, str, list[Union[int, str]]]]]): The index.
            index_2 (Union[int, list[Union[int, str, list[Union[int, str]]]]]): The index.

        Returns:
            float: Correlation function from the samples at the sites of index_1 and index_2.
        """
        return self.__base_solver.calculate_correlation(index_1, index_2)

    def calculate_correlation_list(
        self,
        origin: Union[int, list[Union[int, str, list[Union[int, str]]]]],
        index_list: Union[
            list[int], list[list[Union[int, str, list[Union[int, str]]]]]
        ],
    ) -> list[float]:
        """Calculate correlation function list.

        Args:
            origin (Union[int, list[Union[int, str, list[Union[int, str]]]]]): The index.
            index_list (Union[list[int], list[list[Union[int, str, list[Union[int, str]]]]]]): The index list.

        Returns:
            list[float]: Correlation function list from the samples at the sites of origin and index in index_list.
        """
        return self.__base_solver.calculate_correlation_list(origin, index_list)

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
    def num_threads(self) -> int:
        self.get_num_threads()

    @num_threads.setter
    def num_threads(self, num_threads: int) -> None:
        self.set_num_threads(num_threads=num_threads)

    @property
    def beta(self) -> float:
        return self.get_inverse_temperature()

    @beta.setter
    def beta(self, beta: float) -> None:
        self.set_inverse_temperature(beta=beta)

    @property
    def model(self) -> ModelType:
        return self.__model
