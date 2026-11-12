"""Line of sight from a vehicle in a street canyon to a satellite, and its correlation between two vehicles.

Street along the x axis, walls at y = 0 and y = W, facade heights H drawn per building from a lognormal law
(median H50, log-sd sH), building length Lb along the street. A vehicle antenna at lateral position y and
height h sees a satellite at elevation e and azimuth a (measured from the street axis) if the ray clears the
facade it crosses: the ray reaches the wall on the side sign(sin a) after the horizontal distance
s = d_perp / |sin a|, with d_perp = W - y or y, and clears it iff H < h + s tan e.

    P_LOS(e, a | y) = F_H(h + d_perp(y) tan e / |sin a|),     P_LOS(e, a) = E_y[ P_LOS(e, a | y) ]

The facade that blocks one vehicle is crossed by the ray of a neighbor at x + dx only if the two crossing
points lie on the same building; then the two outcomes share one H and are positively correlated.
"""
import numpy as np
from scipy import stats


class Street:
    def __init__(self, W=20.0, H50=18.0, sH=0.5, Lb=30.0, h=1.5, lanes=(4.0, 7.5, 12.5, 16.0)):
        self.W, self.H50, self.sH, self.Lb, self.h = W, H50, sH, Lb, h
        self.lanes = np.asarray(lanes)

    def F(self, z):
        z = np.asarray(z, float)
        return np.where(z <= 0, 0.0, stats.norm.cdf(np.log(np.maximum(z, 1e-9) / self.H50) / self.sH))

    def dperp(self, y, a):
        return np.where(np.sin(a) >= 0, self.W - y, y)

    def p_los_y(self, e, a, y):
        sa = np.abs(np.sin(a))
        z = self.h + self.dperp(y, a) * np.tan(e) / np.maximum(sa, 1e-9)
        return np.where(sa < 1e-9, 1.0, self.F(z))

    def p_los(self, e, a):
        """Lane-averaged LOS probability (vehicles spread evenly over the lanes)."""
        return float(np.mean([self.p_los_y(e, a, y) for y in self.lanes]))

    def p_both(self, e, a, y1, y2, dx):
        """P(both vehicles in LOS) for two vehicles dx apart along the street at lanes y1, y2."""
        sa, ca = np.abs(np.sin(a)), np.cos(a)
        z1 = self.h + self.dperp(y1, a) * np.tan(e) / max(sa, 1e-9)
        z2 = self.h + self.dperp(y2, a) * np.tan(e) / max(sa, 1e-9)
        # along-street positions where the rays hit the wall
        x1 = self.dperp(y1, a) / max(sa, 1e-9) * ca
        x2 = dx + self.dperp(y2, a) / max(sa, 1e-9) * ca
        gap = abs(x2 - x1)
        same = max(0.0, 1 - gap / self.Lb)       # probability the two crossings fall on one building
        p_same = self.F(min(z1, z2))              # one facade height must clear both rays
        p_diff = self.F(z1) * self.F(z2)
        return same * p_same + (1 - same) * p_diff

    def sample(self, rng, n_buildings):
        return np.exp(np.log(self.H50) + self.sH * rng.standard_normal(n_buildings))


def mc_los(street, e, a, y, x, H, x0=0.0):
    """Exact-geometry LOS for vehicles at (x, y): which building facade the ray crosses, and its height."""
    sa = np.abs(np.sin(a))
    s = street.dperp(y, a) / np.maximum(sa, 1e-9)
    xc = x + s * np.cos(a)
    idx = np.floor((xc - x0) / street.Lb).astype(int) % len(H)
    return H[idx] < street.h + s * np.tan(e)
