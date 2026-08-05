import numpy as np
import warnings
from typing import Callable, Optional

def explicit_fd_python(
    dx: float,
    dt: float,
    M: int,
    Nminus: int,
    Nplus: int,
    pay_off: Callable[[np.ndarray], np.ndarray],
    u_m_inf: Callable[[float, float], float],
    u_p_inf: Callable[[float, float], float],
    values: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Explicit forward-difference solver (1D) for u_t = u_xx (diffusion-like term).
    Params:
      dx, dt: space/time steps
      M: number of time steps
      Nminus, Nplus: integer spatial grid bounds (inclusive)
      pay_off: function pay_off(x_array) returning initial u(x,0)
      u_m_inf: left boundary function u(x_left, tau)
      u_p_inf: right boundary function u(x_right, tau)
      values: optional output array to fill (length Nplus-Nminus+1)
    Returns:
      1D numpy array of solution at final time (size = Nplus - Nminus + 1),
      ordered corresponding to x = np.arange(Nminus, Nplus+1) * dx
    Note:
      Stability check: explicit scheme requires a = dt/dx^2 <= 0.5 for pure diffusion.
    """
    a = dt / (dx * dx)
    if a > 0.5:
        warnings.warn(
            f"Explicit scheme may be unstable: a = dt/dx^2 = {a:.6g} > 0.5. "
            "Consider reducing dt or increasing dx.",
            UserWarning,
        )

    x = np.arange(Nminus, Nplus + 1, dtype=float) * dx
    size = x.size

    # initial condition
    oldu = pay_off(x).astype(float)
    if oldu.shape != (size,):
        raise ValueError("pay_off must return a 1D array of length Nplus-Nminus+1")

    newu = np.empty_like(oldu)

    # time-stepping
    for m in range(1, M + 1):
        tau = m * dt

        # boundary values at current time
        newu[0] = u_m_inf(x[0], tau)
        newu[-1] = u_p_inf(x[-1], tau)

        if size > 2:
            # vectorized interior update for n = Nminus+1 .. Nplus-1
            newu[1:-1] = oldu[1:-1] + a * (oldu[:-2] - 2.0 * oldu[1:-1] + oldu[2:])
        elif size == 2:
            # only boundaries and no interior points: just keep boundaries
            pass
        else:
            # single grid point: update is trivial if needed
            newu[0] = u_m_inf(x[0], tau)  # or keep oldu[0]

        # prepare for next step
        oldu[:] = newu

    if values is None:
        return oldu.copy()
    else:
        if values.shape != oldu.shape:
            raise ValueError("Provided `values` array has wrong shape")
        values[:] = oldu
        return values
