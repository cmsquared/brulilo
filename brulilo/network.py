"""
A Network contains a list of Isotopes, and a list of Reactions that
link the isotopes together.
"""


class Network(object):
    def __init__(self, isotopes, reactions):
        self.isotopes = list(isotopes)
        self.reactions = list(reactions)

        for rxn in self.reactions:
            rxn.update_rxn_vector(isotopes)

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
