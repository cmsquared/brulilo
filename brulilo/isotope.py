"""
An Isotope is a container holding a species mass, atomic number, etc.
"""
import numpy as np

# periodic table, somewhat broken up appropriately
isotope_lut = ['H', 'He',
               'Li', 'Be',
               'B', 'C', 'N', 'O', 'F', 'Ne',
               'Na', 'Mg',
               'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
               'K', 'Ca',
               'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co',  # break
               'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
               'Rb', 'Sr',
               'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh',  # break
               'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
               'Cs', 'Ba',
               # Lanthanide Series
               'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu',  # break
               'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
               'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt',  # break
               'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
               'Fr', 'Ra',
               # Actinide Series
               'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am',  # break
               'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr',
               'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',  # break
               'Rg', 'Cn', 'Uut', 'Fl', 'Uup', 'Lv', 'Uus', 'Uuo']
# this is useful for going from a species name to a Z value
Zdict = {}
for i, species in enumerate(isotope_lut):
    Zdict[species] = i+1


class Isotope(object):
    # some properties for plotting the isotope and it's box
    _width = 0.9
    _label_pad = 0.2
    _box_size = _width - 2*_label_pad

    def __init__(self, name, mass, ebin):
        self.A = mass
        self.Z = Zdict[name.capitalize()]
        self.B = ebin
        self.symbol = isotope_lut[self.Z - 1]
        self._plot_nz = np.array([self.A-self.Z,
                                  self.Z], dtype='int')

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
