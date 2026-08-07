"""EAF slag-metal equilibrium model (Ban-ya quadratic formalism).

Solves the activity coefficients of EAF slag components at a given temperature,
including automatic Fe/Cr redox splitting and slag-metal oxygen-potential balance.
"""

__version__ = "0.1.0"

from .io import InputError, SlagModelConfig, build_report, load_input, parse_input
from .rsm import SlagSplit
from .solver import NonConvergenceError, SolveResult, solve, solve_redox_fixed_po2

__all__ = [
    "InputError",
    "NonConvergenceError",
    "SlagModelConfig",
    "SlagSplit",
    "SolveResult",
    "build_report",
    "load_input",
    "parse_input",
    "solve",
    "solve_redox_fixed_po2",
]
