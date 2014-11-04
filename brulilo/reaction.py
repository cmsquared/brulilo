"""
Storage for a Reaction object, which belongs to a network and acts on
the specified isotopes.  The actual reaction rate is determined by a
reaction function, which is likely a construct for using ReacLib, but
could be anything that takes in a temperature and density.
"""
import numpy as np
import collections
import sys

from isotope import Isotope
from brulilo.util import reaclib as rl


class Reaction(object):
    def __init__(self, species, react_func=None):
        """
        species is a list of tuples, one for each isotope that participates
        in the reaction.  The first element of each tuple in species gives
        the stochiometric coefficient, and the second element is the Isotope.
        """
        coefficients, isotopes = zip(*species)
        self.isotopes = list(isotopes)
        self.coeffs = coefficients
        self.rate = react_func

    def __str__(self):
        lhs = []
        rhs = []
        for coeff, isotope in zip(self.coeffs, self.isotopes):
            if coeff < 0:
                lhs.extend([str(isotope), ]*(-coeff))
            else:
                rhs.extend([str(isotope), ]*coeff)
        lhs = ' + '.join(lhs)
        rhs = ' + '.join(rhs)
        return lhs + ' --> ' + rhs

    def update_rxn_vector(self, network_isotopes):
        """
        Each reaction has a 'direction' in a multidimensional space
        defined by the total isotopes in a network.  We use this
        to determine how a particular reaction acts within a network.
        """
        # sort the isotopes in some predictable fashion
        network_isotopes.sort()
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

    def plot_on(self, fig):
        """
        Plop the reaction onto a figure.
        """
        pass


class ReacLibReaction(Reaction):
    """
    For ReacLib[1] reaction rate objects.  These will be stored in an XML
    file that we only want to parse once, so we pass just the name of the
    isotope -- NOT the Isotope object -- to the baseclass.

    [1] https://groups.nscl.msu.edu/jina/reaclib/db/
    """
    def __init__(self, nucRxnString):
        """
        nucRxnString is a reaction descriptor in the form of a(b,c)d
        """
        self.nucRxnString = nucRxnString
        # build target, reactants and products
        self._parse_rate_string()

        # separate into unique values
        isotopes_cnt = collections.Counter(self.reactants + self.products)
        uniq_reactants = set(self.reactants)
        uniq_products = set(self.products)
        species = [(-1*isotopes_cnt[isotope], isotope)
                   for isotope in uniq_reactants]
        species += [(isotopes_cnt[isotope], isotope)
                    for isotope in uniq_products]

        super(ReacLibReaction, self).__init__(species)

    def _parse_rate_string(self):
        """
        Takes a nucRxnString in the form of He4(aa,g)C12 and parses it into
        the target nucleus - He4 -, the reactants - [He4,aa] -, and the
        products - [g,C12] - useful for searching the ReacLib database.
        Handling of special characters (e.g. aa, g) is done via the
        sanitize_species method from the util.reaclib module.
        """
        # tokenize the nucRateString and apply any sanitizations
        target = rl.sanitize_species(self.nucRxnString.split('(')[0])
        internals = self.nucRxnString.split('(')[1].split(')')[0].split(',')
        reactants, products = [rl.sanitize_species(spec) for spec in internals]
        endState = rl.sanitize_species(self.nucRxnString.split(')')[1])

        # this is of the form a + b -> c + d, with special character handling
        rateString = rl.form_rate_string(target, reactants, products, endState)

        # split this into reactants and products, for easier searching of the
        # various permutations in the ReacLib database
        self.reactants = rateString.split(rl.TO)[0].split(rl.PLUS)
        self.products = rateString.split(rl.TO)[1].split(rl.PLUS)

        self.target = target

    def build_rxn_rate(self, xml_root):
        """
        xml_root is the root node of the XML document containing all the
        reaction rate information in the ReacLib database.
        """
        # xpath_str = '//reaction[%s]'
        # reactants = ['reactant="%s"' % reactant.lower()
        #              for reactant in self.reactants]
        # products = ['product="%s"' % product.lower()
        #             for product in self.products]
        # xpath_str = xpath_str % ' and '.join(reactants + products)
        # search for reactions that have these specific reactants and products
        possible_reactions = xml_root.xpath(_build_xpath(self.reactants,
                                                         self.products))
        # if possible_reactions is an empty list, then this might be
        # a reverse rate, which isn't stored in the database.  need to
        # construct the reverse rate, if possible
        if not possible_reactions:
            possible_reactions = xml_root.path(_build_xpath(self.products,
                                                            self.reactants))
            if not possible_reactions:
                errStr = "%s\n or it's reverse not found in ReacLib database."
                raise RuntimeError(errStr % str(self))

            # FIXME - implement this
            _build_reverse()

        # possible_reactions is NOT a list of reactions with ONLY these
        # products and reactants, so we need to filter some out by making
        # sure the NUMBER of products and reactants match, respectively
        rxns_same_no_reactants = [rxn for rxn in possible_reactions if
                                  int(rxn.xpath('count(.//reactant)')) ==
                                  len(self.reactants)]
        this_reaction = [rxn for rxn in rxns_same_no_reactants if
                         int(rxn.xpath('count(.//product)')) ==
                         len(self.products)]
        if len(this_reaction) != 1:
            errStr = "Found multiple reactions for \n%s\n" % str(self)
            errStr += xpath_str
            raise RuntimeError(errStr)
        this_reaction = this_reaction[0]
        # now read in the recommended fit parameters
        aFactors = {}
        for a_param in ['a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7']:
            aFactors[a_param] = [float(a_val.text) for a_val in
                                 this_reaction.xpath('//%s' % a_param)]

        self.rate = rl.build_rate_function(aFactors)


def _build_xpath(reactants, products):
    # use XPath syntax to search for this specific rate based on
    # the reactants and products
    xpath_str = '//reaction[%s]'
    reacts = ['reactant="%s"' % reactant.lower()
              for reactant in reactants]
    prods = ['product="%s"' % product.lower()
             for product in products]
    return xpath_str % ' and '.join(reacts + prods)

def _build_reverse():
    """
    Need to implement reverse rates here
    """
    pass
