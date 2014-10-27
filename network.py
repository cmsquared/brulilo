"""
A Network contains a list of Isotopes, and a list of Reactions that
link the isotopes together.
"""


class Network(object):
    def __init__(self, isotopes, reactions):
        self.isotopes = list[isotopes]
        self.reactions = list[reactions]

        for rxn in self.reactions:
            rxn.update_rxn_vector(isotopes)
