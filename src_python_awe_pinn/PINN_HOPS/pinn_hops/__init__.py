"""PINN solver for the 2D two-layer grating problem of Kehoe & Nicholls (J. Sci. Comput. 2024), eqs. (6),
and comparison tools against the HOPS/AWE solver (HOPS_Python).

GratingPINN / FieldNet (the PyTorch PINN) are imported lazily, so the numpy-only parts (problem, lsq_pinn,
hops_reference) -- all that HOPS_PINN_Hybrid and HOPS_AWE_PINN use -- work without torch installed."""
from .problem import Grating2D, outgoing_sqrt, PROFILES


def __getattr__(name):
    if name in ('GratingPINN', 'FieldNet'):
        from . import pinn
        return getattr(pinn, name)
    raise AttributeError(f"module 'pinn_hops' has no attribute {name!r}")
