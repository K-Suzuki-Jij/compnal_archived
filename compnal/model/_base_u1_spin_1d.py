import base_compnal
from functools import singledispatch


class U1Spin1D(base_compnal.model._U1Spin_1D):
    """One-dimensional spin systems with the U(1) symmetry.
    
    Args:
        system_size (int): The system size.
        spin (float): Magnitude of spin.

    
    """
    def __init__(self, system_size: int = 0, spin: float = 0.5, total_sz = 0) -> None:
        super().__init__(system_size, spin, total_sz)

    def add_onsite_potential(self, m, site: int, value: float=1.0):
        super()._add_onsite_potential(m, site, value)

    def add_interaction(self, *args, **kwargs):
        super()._add_interaction(*args, **kwargs)

    def print_basis_onsite(self) -> str:
        magnitude_spin = super()._get_magnitude_spin()
        dim_onsie = super()._get_dim_onsite()
        for row in range(dim_onsie):
            print("row ", row, ": |Sz={:+}".format(magnitude_spin - row), ">", sep='')

    @singledispatch
    def calculate_target_dim(self, *args, **kwargs):
        return super()._calculate_target_dim(*args, **kwargs)

    @staticmethod
    @calculate_target_dim.register
    def _calculate_target_dim(system_size: int, spin: float, total_sz: float) -> int:
        if system_size < 0:
            raise RuntimeError("system_size must be non-negative integer")

        if not base_compnal.model._U1Spin_1D.is_valid_q_number(system_size, spin, total_sz):
            return 0
        
        magnitude_2spin = int(2*spin)
        total_2sz = int(2*total_sz)
        max_total_2sz = system_size*magnitude_2spin

        dim = [[0]*(max_total_2sz + 1) for _ in range(system_size)]

        for s in range(-magnitude_2spin, magnitude_2spin + 1, 2):
            dim[0][(s + magnitude_2spin)//2] = 1

        for site in range(1, system_size, 1):
            for s in range(-magnitude_2spin, magnitude_2spin + 1, 2):
                for s_prev in range(-magnitude_2spin*site, magnitude_2spin*site + 1, 2):
                    a = dim[site    ][(s + s_prev + magnitude_2spin*(site + 1))//2]
                    b = dim[site - 1][(s_prev + magnitude_2spin*site)//2]
                    dim[site][(s + s_prev + magnitude_2spin*(site + 1))//2] = a + b

        return dim[system_size - 1][(total_2sz + max_total_2sz)//2]

    @property
    def dim_onsite(self):
        return super()._get_dim_onsite()

    @property
    def total_sz(self):
        return super()._get_total_sz()

    @total_sz.setter
    def total_sz(self, total_sz: float()):
        super()._set_total_sz(total_sz)
    
    @property
    def spin(self):
        return super()._get_magnitude_spin()

    @spin.setter
    def spin(self, spin: float):
        super()._set_magnitude_spin(spin)

    @property
    def system_size(self):
        return super()._get_system_size()

    @system_size.setter
    def system_size(self, system_size):
        super()._set_system_size(system_size)

    @property
    def sx(self):
        return super()._get_onsite_operator_sx()

    @property
    def isy(self):
        return super()._get_onsite_operator_isy()

    @property
    def sz(self):
        return super()._get_onsite_operator_sz()

    @property
    def sp(self):
        return super()._get_onsite_operator_sp()

    @property
    def sm(self):
        return super()._get_onsite_operator_sm()

    @staticmethod
    def make_onsite_operator_sx(spin: float):
        return base_compnal.model._U1Spin_1D._make_onsite_operator_sx(spin)

    @staticmethod
    def make_onsite_operator_isy(spin: float):
        return base_compnal.model._U1Spin_1D._make_onsite_operator_isy(spin)

    @staticmethod
    def make_onsite_operator_sz(spin: float):
        return base_compnal.model._U1Spin_1D._make_onsite_operator_sz(spin)

    @staticmethod
    def make_onsite_operator_sp(spin: float):
        return base_compnal.model._U1Spin_1D._make_onsite_operator_sp(spin)

    @staticmethod
    def make_onsite_operator_sm(spin: float):
        return base_compnal.model._U1Spin_1D._make_onsite_operator_sm(spin)
