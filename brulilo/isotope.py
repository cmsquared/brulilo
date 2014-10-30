"""
An Isotope is a container holding a species mass, atomic number, etc.
"""
import numpy as np
from util.reaclib import get_Z_A, isotope_lut, Zdict

isotope_registry = {}


class Isotope(object):
    # some properties for plotting the isotope and it's box
    _width = 0.9
    _label_pad = 0.2
    _box_size = _width - 2*_label_pad

    def __new__(cls, name, mass=None, **kwargs):
        if mass is None:
            Z, A = get_Z_A(name)
        else:
            A = mass
            Z = Zdict[name.capitalize()]
        registered_name = "%s%d" % (isotope_lut[Z-1], A)
        if registered_name in isotope_registry:
            return isotope_registry[registered_name]
        return super(Isotope, cls).__new__(cls, name, mass, **kwargs)

    def __init__(self, name, mass=None, ebin=None):
        """
        You can specify the name of the isotope (e.g. "He") and the mass
        (e.g. 4) separately, or as a single string (e.g. "He4").
        """
        if mass is None:
            self.Z, self.A = get_Z_A(name)
        else:
            self.A = mass
            self.Z = Zdict[name.capitalize()]
        self.B = ebin
        # special neutron case
        if self.Z == 0:
            self.symbol = 'n'
        else:
            self.symbol = isotope_lut[self.Z - 1]
        self._plot_nz = np.array([self.A-self.Z,
                                  self.Z], dtype='int')

        # probably a more meta way of doing this, but register this class
        isotope_registry[str(self)] = self

    def __hash__(self):
        # just so we can make a set...
        return isotope_registry.keys().index(str(self))

    def __cmp__(self, other):
        return cmp(hash(self), hash(other))

    def __str__(self):
        return "%s%d" % (self.symbol, self.A)

    def _plot_build_label(self):
        """
        Provide a nice LaTeX-formatted label.  Useful for plotting.
        """
        baseStr = r'$^{%d}\mathrm{%s}$'
        return baseStr % (self.A, self.symbol)

    def plot_label(self, fig):
        """
        Plop the label onto a figure.
        """
        fig.gca().text(self._plot_nz[0] + 0.5*self._width,
                       self._plot_nz[1] + 0.5*self._width,
                       self._plot_build_label(),
                       color='black',
                       # family='sans-serif',
                       fontsize=14,
                       ha='center', va='center',
                       # make this one of the last things drawn
                       zorder=100)

    def plot_patch(self, fig):
        """
        Plop the box for this Isotope on a figure.
        """
        from matplotlib.patches import FancyBboxPatch
        fig.gca().add_patch(FancyBboxPatch(self._plot_nz + self._label_pad,
                                           self._box_size, self._box_size,
                                           ec='black', fc='none',
                                           boxstyle=('round,pad=%s' %
                                                     self._label_pad),
                                           ls='solid'))
