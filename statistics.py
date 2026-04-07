import numpy as np


def mandel_paule(values: list, uncertainties: list) -> tuple:
    """Compute the Mandel-Paule weighted mean and its uncertainty.
    References:
        Paule, R.C. and Mandel, J., "Consensus Values and Weighting
        Factors", J. Res. Natl. Bur. Stand., 87(5), 377-385, 1982.

    Parameters:
        values: Measurement values x_i.
        uncertainties: Standard uncertainties u_i (must be > 0).

    Returns:
        (weighted_mean, u_weighted_mean, dark_uncertainty)

    Example:
        >>> vals = [10.1, 10.4, 9.8, 10.3, 10.0]
        >>> unc = [0.10, 0.20, 0.15, 0.12, 0.18]
        >>> mean, u_mean, du = mandel_paule(vals, unc)
    """
    tol = 1e-8
    x = np.asarray(values, dtype=float)
    u2 = np.asarray(uncertainties, dtype=float) ** 2
    nu = len(x) - 1

    def calc(dark2):
        wgt = 1.0 / (u2 + dark2)
        s_w = wgt.sum()
        x_bar = (wgt * x).sum() / s_w
        return x_bar, wgt, s_w, (wgt * (x - x_bar) ** 2).sum()

    xbar, w, sw, chi2 = calc(0.0)

    du2 = 0.0
    if chi2 > nu:
        lo, hi = 0.0, float(max(np.var(x, ddof=1), u2.max()) * 10.0)
        while calc(hi)[3] >= nu:
            hi *= 10.0
        for i in range(100):
            mid = 0.5 * (lo + hi)
            if calc(mid)[3] > nu:
                lo = mid
            else:
                hi = mid
            if (hi - lo) < tol * max(mid, 1e-30):
                # print(f"Iterations: {i+1}")
                break

        du2 = 0.5 * (lo + hi)
        xbar, w, sw, _ = calc(du2)

    return xbar, 1.0 / np.sqrt(sw), np.sqrt(du2)
