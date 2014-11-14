"""
An Isotope is a container holding a species mass, atomic number, etc.
"""
import numpy as np
from scipy.interpolate import interp1d

from util.species import get_Z_A, element_lut, Zdict
from util.constants import MeV2erg
from util.webnucleo import webnucleo

isotope_registry = {}


class Isotope(object):
    # some properties for plotting the isotope and it's box
    _width = 0.9
    _label_pad = 0.2
    _box_size = _width - 2*_label_pad

    def __new__(cls, name, **kwargs):
        Z, A = get_Z_A(name)
        if Z != 0:
            registered_name = "%s%d" % (element_lut[Z], A)
        else:  # neutron
            registered_name = "n"
        if registered_name not in isotope_registry:
            this_isotope = super(Isotope, cls).__new__(cls, name, **kwargs)
            isotope_registry[registered_name] = this_isotope

        return isotope_registry[registered_name]

    def __init__(self, name, pbar=None):
        """
        You specify the name of the isotope and the mass as a single string
        (e.g. "He4").  Nuclear properties, like mass_excess, will be looked up
        in a nuclear data file from ReacLib/Webnucleo.

        pbar can be a progressbar instance, useful for ticking along for large
        networks to show progress.
        """
        self.Z, self.A = get_Z_A(name)
        self.symbol = element_lut[self.Z]
        self._plot_nz = np.array([self.A-self.Z,
                                  self.Z], dtype='int')

        # look up the nuclear data and build the partition function
        if pbar is not None:
            pbar.update(str(self))
        self._build_nuclear_data()

    def _build_nuclear_data(self):
        # find my entry in the Webnucleo data file
        my_data = webnucleo.get_isotope_data(self)
        # mass excess and spin
        # TODO -- make this work with nuclides with different states like Al26
        for attr in ["mass_excess", "spin"]:
            setattr(self, attr, float(my_data.find(attr).text))
        # binding energy
        self.binding_energy = (self.Z*webnucleo.proton_mass_excess +
                               (self.A-self.Z)*webnucleo.neutron_mass_excess -
                               self.mass_excess)
        # swap to cgs units
        self.binding_energy *= MeV2erg
        self.mass_excess *= MeV2erg
        # get the partition table entry
        ptable = my_data.find("partf_table")
        if ptable is None:
            # partition table data does not exist, so the partition function
            # is just the number of states of the ground state
            self.partition_function = lambda temperature: 2*self.spin + 1
        else:
            # now the actual data: (t9, log10(f)) pairs, where f is related
            # to the partition function -- see Webnucleo documentation
            part_t9 = map(float, [t9.text for t9 in ptable.xpath("point/t9")])
            part_lf = map(float, [lf.text
                                  for lf in ptable.xpath("point/log10_partf")])
            self.partition_function = self._build_partition_function(part_t9,
                                                                     part_lf)

    def _build_partition_function(self, table_t9, table_logf):
        fit = interp1d(table_t9, table_logf, kind='cubic')
        minT9 = min(table_t9)
        maxT9 = max(table_t9)

        def _part_function(isotope, temperature):
            t9 = temperature / 1e9
            # no extrapolation
            t9 = min(maxt9, max(mint9, t9))
            return (2*self.spin + 1) * 10**fit(t9)

        return _part_function

    def __hash__(self):
        # just so we can make a set...
        return isotope_registry.keys().index(str(self))

    def __cmp__(self, other):
        # just so we can make a set...
        return cmp(hash(self), hash(other))

    def __str__(self):
        my_str = self.symbol
        if my_str != 'n':
            my_str += str(self.A)
        return my_str

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
