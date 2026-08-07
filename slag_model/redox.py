"""Redox splitting: Fe (Ban-ya Eq. 24) and Cr (exact cubic)."""

from __future__ import annotations

import math

from .data import (
    BANYA_CONSTANT,
    BANYA_PO2_COEF,
    BANYA_T_COEF,
)


def fe_ratio_ban_ya(T: float, P_O2: float, gamma_FeO: float, gamma_FeO1_5: float) -> float:
    """Fe3+/Fe2+ ratio from Ban-ya Eq. 24 at fixed P_O2 (PLAN.md 6.4d, 7.5).

    log10(Fe3+/Fe2+) = 6625/T - 2.77 + 0.25*log10(P_O2)
                        + log10(gamma_FeO) - log10(gamma_FeO1.5)
    """
    log10_r = (
        BANYA_T_COEF / T
        + BANYA_CONSTANT
        + BANYA_PO2_COEF * math.log10(P_O2)
        + math.log10(gamma_FeO)
        - math.log10(gamma_FeO1_5)
    )
    return 10.0**log10_r


def cr_ratio(rhs: float) -> float:
    """Solve r^3/(1+r) = rhs for r > 0 by bisection (PLAN.md 7.6).

    f(r) = r^3/(1+r) - rhs is strictly increasing for r > 0, so the root in
    [1e-6, 1e6] is unique. Bisection to machine precision in 100 iterations.
    """
    if rhs <= 0.0:
        return 0.0
    lo = 1e-6
    hi = 1e6
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if mid * mid * mid / (1.0 + mid) < rhs:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)
