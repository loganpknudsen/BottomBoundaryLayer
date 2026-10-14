"""Floquet stability analysis of PSI in the bottom boundary layer (manuscript Section 2).

Solves Eq. 2.26 with the coefficients of Eq. 2.27, derived from Eq. 2.25 without the further
small-angle reduction of Eq. 2.28. Note that the A0 printed as Eq. 2.27a in the May 2026
preprint is missing a +S^2 gamma term inside the bracket; the form used here is

    A0 = -[mu0 (S^2/theta - theta) + S^2 (1 - gamma)/theta^2 - S^2 gamma + theta^2] / (1 + S^2)
    A1 = mu0 + (S^2 - theta^2) / (theta (1 + S^2))

For this system the linear TKE budget (Eq. 2.31) is exact, so growth rates, perturbation slopes
and budget terms all come from one set of equations. Time is t_hat = f* t with
f* = f sqrt(1 + S^2), so growth rates are in units of f*.
"""
import numpy as np
from scipy.integrate import simpson, solve_ivp
from scipy.optimize import minimize_scalar

THETA = 0.1
F = 1e-4
PHI = np.pi / 2           # wave phase, as in the simulations' (H - z) convention
TAU = 2 * np.pi           # one wave period in t_hat
RTOL, ATOL = 1e-9, 1e-12


class Flow:
    """Background parameters for one (S, gamma, delta)."""

    def __init__(self, S, gamma, delta, theta=THETA, f=F):
        self.S, self.gamma, self.delta, self.theta, self.f = S, gamma, delta, theta, f
        self.S2 = S * S
        self.N2 = self.S2 * f**2 / theta**2
        self.fstar = f * np.sqrt(1 + self.S2)
        self.Lam = self.N2 * theta * gamma / f          # geostrophic shear Lambda
        self.amp = delta * self.Lam / self.fstar        # wave amplitude of mu(t)

    def A(self, mu0):
        S2, th, g = self.S2, self.theta, self.gamma
        A0 = -(mu0 * (S2 / th - th) + S2 * (1 - g) / th**2 - S2 * g + th**2) / (1 + S2)
        A1 = mu0 + (S2 - th**2) / (th * (1 + S2))
        return A0, A1

    def mu(self, t, mu0):
        return mu0 + self.amp * np.sin(t + PHI)


def monodromy(flow, mu0s):
    """Monodromy matrices over one period for an array of mu0. Returns shape (n, 2, 2),
    columns being the solutions from (Psi, zeta) = (1, 0) and (0, 1)."""
    mu0s = np.atleast_1d(np.asarray(mu0s, dtype=float))
    A0, A1 = flow.A(mu0s)

    def rhs(t, y):
        y = y.reshape(4, -1)
        mu = flow.mu(t, mu0s)
        a = A0 - A1 * mu
        r = 1.0 / (1.0 + mu**2)
        return np.concatenate([a * y[1], r * y[0], a * y[3], r * y[2]])

    n = mu0s.size
    y0 = np.concatenate([np.ones(n), np.zeros(n), np.zeros(n), np.ones(n)])
    yT = solve_ivp(rhs, (0.0, TAU), y0, method="DOP853", rtol=RTOL, atol=ATOL).y[:, -1].reshape(4, n)
    M = np.empty((n, 2, 2))
    M[:, 0, 0], M[:, 1, 0], M[:, 0, 1], M[:, 1, 1] = yT[0], yT[1], yT[2], yT[3]
    return M


def exponents(M):
    """Floquet exponent with the largest real part for each monodromy matrix:
    returns (sigma, |Im sigma|) in units of f*."""
    lam = np.linalg.eigvals(M).astype(complex)
    ex = np.log(lam) / TAU
    i = np.argmax(ex.real, axis=-1)
    best = np.take_along_axis(ex, i[..., None], axis=-1)[..., 0]
    return best.real, np.abs(best.imag)


def fastest_mode(flow, m_grid):
    """Maximise the growth rate over mu0: a scan of m_grid, then a bounded refinement between
    the neighbours of the best grid point. Returns (sigma, frequency, mu0)."""
    sig, _ = exponents(monodromy(flow, m_grid))
    i = int(np.argmax(sig))
    lo, hi = m_grid[max(i - 1, 0)], m_grid[min(i + 1, len(m_grid) - 1)]
    res = minimize_scalar(lambda m: -exponents(monodromy(flow, m))[0][0],
                          bounds=(lo, hi), method="bounded", options={"xatol": 1e-5})
    m = res.x if -res.fun >= sig[i] else m_grid[i]
    s, fr = exponents(monodromy(flow, m))
    return s[0], fr[0], m


BUDGET_TERMS = ("GSP", "ASPx", "ASPy", "BPx", "BPz")


def budget(flow, mu0, ic="floquet", nt=2001):
    """Eq. 2.31 terms integrated over one period, divided by the change in k_hat.

    ic="floquet" starts from the eigenvector of the growing multiplier, so the result does not
    depend on the starting phase; ic="ic01" starts from (Psi, zeta) = (0, 1) as in the preprint.
    Returns a dict of the five fractions plus 'closure' = (sum of terms)/dk - 1 and 'dk'.
    """
    nan = {k: np.nan for k in BUDGET_TERMS + ("closure", "dk")}
    if not (0 < flow.gamma < 1):
        return nan
    if ic == "floquet":
        lam, vec = np.linalg.eig(monodromy(flow, mu0)[0])
        j = int(np.argmax(np.abs(lam)))
        if np.iscomplexobj(lam) and lam[j].imag != 0:   # complex pair: no growing real mode
            return nan
        y0 = np.real(vec[:, j])
    else:
        y0 = np.array([0.0, 1.0])
    A0, A1 = flow.A(mu0)
    rhs = lambda t, y: [(A0 - A1 * flow.mu(t, mu0)) * y[1], y[0] / (1 + flow.mu(t, mu0) ** 2)]
    t = np.linspace(0.0, TAU, nt)
    Psi, zeta = solve_ivp(rhs, (0.0, TAU), y0, t_eval=t, method="DOP853", rtol=RTOL, atol=ATOL).y

    f, th, g, S2, N2, Lam = flow.f, flow.theta, flow.gamma, flow.S2, flow.N2, flow.Lam
    mu = flow.mu(t, mu0)
    n = flow.fstar / Lam
    w = Psi / (1 + mu**2)
    u = mu * w * n
    v = -f * (g * S2 / th - th + mu0) * zeta / Lam
    b = -N2 * th * (mu0 + (1 - g) / th) * zeta / (N2 * (1 - g))
    k = 0.5 * (u**2 + v**2 + n**2 * w**2)
    Ri = N2 * (1 - g) / Lam**2
    d = flow.delta
    terms = {
        "GSP": -v * w,
        "ASPx": -d * np.cos(t + PHI) * u * w,
        "ASPy": d * f / flow.fstar * np.sin(t + PHI) * v * w,
        "BPx": (1 - g) / (g * np.sqrt(1 + S2)) * u * b,
        "BPz": Ri * w * b,
    }
    dk = k[-1] - k[0]
    out = {name: simpson(val, x=t) / dk for name, val in terms.items()}
    out["closure"] = sum(out[nm] for nm in BUDGET_TERMS) - 1
    out["dk"] = dk
    return out
