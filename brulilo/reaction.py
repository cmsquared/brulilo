"""
Storage for a Reaction object, which belongs to a network and acts on
the specified isotopes.  The actual reaction rate is determined by a
reaction function, which is likely a construct for using ReacLib, but
could be anything that takes in a temperature and density.
"""
import numpy as np


class Reaction(object):
    def __init__(self, species, react_func):
        """
        species is a list of tuples, one for each isotope that participates
        in the reaction.  The first element of each tuple in species gives
        the stochiometric coefficient, and the second element is the Isotope.
        """
        coefficients, isotopes = zip(*species)
        self.isotopes = list(isotopes)
        self.coeffs = coefficients
        self.func = react_func

    def update_rxn_vector(self, network_isotopes):
        """
        Each reaction has a 'direction' in a multidimensional space
        defined by the total isotopes in a network.  We use this
        to determine how a particular reaction acts within a network.
        """
        vec = np.zeros(len(network_isotopes), dtype='int')
        indices = [i for i,isotope in enumerate(network_isotopes)
                   if isotope in self.isotopes]
        vec[indices] = 1
        for i, isotope in enumerate(self.isotopes):
            try:
                i_net = network_isotopes.index(isotope)
            except IndexError:
                raise RuntimeError("%s not in network" % isotope)
            else:
                vec[i_net] *= self.coeffs[i]
        self.rxn_vector = vec[:]
