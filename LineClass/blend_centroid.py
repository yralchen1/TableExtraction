"""The uncertainty a blend's centroid owes to its calculated intensities.

A line classified to several transitions is fitted as the centroid of its
components: the observed wavenumber is compared with

    x_cg = sum_i w_i * x_i,     w_i = I_i / sum_j I_j,

x_i the wavenumber component i is predicted at and I_i its calculated
intensity.  The w_i are only as good as the I_i, and nothing else in the
pipeline carries their error: the line's own uncertainty is that of the
measurement, and LOPT's centroid model takes the weights as exact.  A blend
whose components lie far apart then pins its levels harder than the data
allow.

If each ln I_i is uncertain by `u_ln`, independently, then since
d x_cg / d ln I_i = w_i * (x_i - x_cg),

    u_I = u_ln * sqrt( sum_i (w_i * (x_i - x_cg))^2 ),

which is added to the line's uncertainty in quadrature wherever the blend
is weighted in a level fit (classify_lines.calc_weights,
make_LOPT_input.write_lines_file).  Two equal components a distance delta
apart give u_I = u_ln * delta / sqrt(8); a lopsided pair gives less, since
the weak component can move the centroid only a little.

This is not the blending factor k(n) of level_positions.py, which measures
the spread of the components about the centroid and which LOPT's centroid
model removes.  u_I is the uncertainty of the centroid's own position.

`u_ln` is `[blends] u_ln_intensity` of the configuration.  0.33 (2026-10-09)
is the upper end of the 68 % interval measured from the leave-one-out
residuals of the 515 blends of iter_hfs that had no registry inflation:
maximum likelihood gave 0.24 (0.14-0.33), the same for plain lines and for
the others.  The upper end was taken because only blends whose centroid
fitted were accepted, so the accepted ones understate the scatter.
"""
import math


def intensity_unc(intensities, positions, u_ln):
    """u_I of a blend whose components have the calculated `intensities` and
    sit at `positions` (cm^-1, or any offsets from a common origin), each
    ln-intensity uncertain by `u_ln`.  0 for fewer than two components, for
    `u_ln` <= 0, or when the intensities add up to nothing."""
    total = sum(intensities)
    if u_ln <= 0 or len(intensities) < 2 or total <= 0:
        return 0.0
    w = [i / total for i in intensities]
    cg = sum(wi * x for wi, x in zip(w, positions))
    return u_ln * math.sqrt(sum((wi * (x - cg)) ** 2
                                for wi, x in zip(w, positions)))
