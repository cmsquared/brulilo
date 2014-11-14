"""
A Network contains a list of Isotopes, and a list of Reactions that
link the isotopes together.
"""
import types
import lxml.etree as etree
import os.path

from reaction import Reaction
from isotope import Isotope
from util.progressbar import IntProgressBar
# import brulilo.util.reaclib as rl
import brulilo


class Network(object):
    def __init__(self, isotopes, reactions):
        self.isotopes = list(isotopes)
        self.reactions = list(reactions)

        # if the isotopes are just the string names, as is the case
        # for from_rxn_file-generated networks, we need to do some
        # updating
        if isinstance(self.isotopes[0], types.StringTypes):
            self._update_with_Isotopes()
            rl.get_nuclear_data(self)
            rl.get_rate_data(self)
#            self._build_rxn_rates()

        for rxn in self.reactions:
            rxn.update_rxn_vector(self.isotopes)

    @classmethod
    def from_rxn_file(cls, rxn_file):
        reactions = []
        with open(rxn_file, 'r') as f:
            nrxns = sum(1 for line in f)
        pbar = IntProgressBar('Building rxns and Isotopes', nrxns)
        with open(rxn_file, 'r') as f:
            for rxn in f:
                reactions.append(Reaction(rxn.strip(), pbar=pbar))
        for reaction in reactions:
            print str(reaction)
            print reaction.isotopes
            print reaction.reactants
            print reaction.products
        isotopes = []
        for reaction in reactions:
            isotopes.extend(reaction.isotopes)

        isotopes = set(isotopes)
        return cls(isotopes, reactions)

    def _update_with_Isotopes(self):
        """
        Turn the strings stored in self.isotopes into actual Isotope objects.
        Furthermore, update each ReacLibReaction to using the Isotopes.
        """
        real_isotopes = []
        for isotope in self.isotopes:
            real_isotopes.append(Isotope(isotope))
        # let's get each Isotope's mass_excess from the nuclear data file
        # this takes a bit
        # fn = os.path.join(os.path.dirname(brulilo.__file__), nuc_data_file)
        # nuc_data_root = etree.parse(fn)
        # for isotope in real_isotopes:
        #     isotope.update_mass_excess(nuc_data_root)
        self.isotopes = real_isotopes

        net_isotopes = {str(isotope): isotope for isotope in self.isotopes}
        for reaction in self.reactions:
            reaction.isotopes = [net_isotopes[isotope] for
                                 isotope in reaction.isotopes]

    def _build_rxn_rates(self):
        """
        Read the reaction fit parameters from the ReacLib database stored
        in the XML file, rxn_data_file; may take a while
        """
        fn = os.path.join(os.path.dirname(brulilo.__file__), rxn_data_file)
        rxn_data_root = etree.parse(fn)
        for reaction in self.reactions:
            reaction.build_rxn_rate(rxn_data_root)

    def pprint(self):
        print 'Isotopes:'
        for isotope in self.isotopes:
            print isotope
        print 'Reactions:'
        for reaction in self.reactions:
            print reaction

    def plot(self, draw_rxns=False):
        """
        Draw a network diagram of the Isotopes and potentially the Reactions.
        """
        import matplotlib.pyplot as plt

        fig = plt.figure()

        for isotope in self.isotopes:
            isotope.plot_label(fig)
            isotope.plot_patch(fig)
        if draw_rxns:
            # do Reaction drawing stuff here
            pass

        # fix the ticks
        ax = fig.gca()
        ax.set_aspect('equal')
        for side in ['right', 'top']:
            ax.spines[side].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.set_xticks([])
        ax.set_yticks([])

        # fix the limits of the plot
        nlim, zlim = self._determine_plot_extent()
        ax.set_xlim(*nlim)
        ax.set_ylim(*zlim)

        # add some labeling
        ax.set_xlabel("$N\ \longrightarrow$", fontsize=16)
        ax.set_ylabel("$Z\ \longrightarrow$", fontsize=16)

        plt.show()

    def _determine_plot_extent(self):
        ns = [isotope.A - isotope.Z for isotope in self.isotopes]
        zs = [isotope.Z for isotope in self.isotopes]
        # use min/max, but with padding
        general_pad = self.isotopes[0]._label_pad
        box_width_pad = self.isotopes[0]._width
        nextent = [min(ns) - general_pad,
                   max(ns) + general_pad + box_width_pad]
        zextent = [min(zs) - general_pad,
                   max(zs) + general_pad + box_width_pad]
        return nextent, zextent
