"""
Storage for a Reaction object, which belongs to a network and acts on
the specified Isotopes.  The actual reaction rate is determined by a
reaction function, which is likely a construct for using ReacLib, but
could be anything that takes in a temperature and density.
"""
import numpy as np
# import collections
# import sys
import re

from util.webnucleo import webnucleo
from util.species import \
    sanitize_species, form_rate_string, PLUS, isotope_lut, leptons
from isotope import Isotope
from util.constants import electron_mass, light_speed

rate_breaker = re.compile("(.*)\((.*),(.*)\)(.*)")

class Reaction(object):

    # some properties that modify some values
    is_reverse = False
    is_weak = False
    is_betaplus = False
    is_electron_capture = False

    def __init__(self, rxnString, pbar=None):
        """
        rxnString is the reaction rate in typical astrophysical notation, 
        including leptons using syntax defined in the README.  For example,
        the reaction

                 p + p --> d + e^+ + nu

        would have a rxnString="p(p,e+nu_e)d" or one could write out the
        nuclides explicitly as rxnString="H1(H1,e+nu_e)H2".

        pbar is a progress bar object, which is useful for tracking the 
        progress of Network creation on large numbers of Reactions.
        """
        self.rxnString = rxnString
        # parse into reactants and products
        self._parse_rate()

        # some reaction rate qualifiers
        self.is_weak = any([species in self.reactants + self.products
                            for species in leptons])
        self.is_betaplus = all([species in self.products
                                for species in ["electron", "neutrino_e"]])
        self.is_electron_capture = (any([species in self.reactants
                                         for species in ["positron",
                                                         "anti-neutrino_e"]])
                                    and
                                    any([species in self.reactants
                                         for species in ["positron",
                                                         "anti-neutrino_e"]]))

        # make the Isotope objects for this reaction
        self.isotope_reactants = [Isotope(reactant) for reactant in
                                  self.reactants
                                  if reactant in isotope_lut]
        self.isotope_products = [Isotope(product) for product in
                                 self.products
                                 if product in isotope_lut]

        # let's find our reaction data in the data file
        if pbar is not None:
            pbar.update(self.rxnString)
        self._build_rxn_data()
        self._build_qvalue()

    def __str__(self):
        # kludge to use the form_rate_string
        return form_rate_string([self.reactants[0]], self.reactants[1:],
                                self.products[1:], [self.products[0]])

    def _parse_rate(self):
        """
        Takes a rxnString in the form of 
        
        He4(aa,g)C12

        and parses it into the target nucleus - He4 -, the reactants -
        [He4,He4,He4] -, and the products - [gamma,C12] - useful for
        searching the ReacLib database.  Handling of special
        characters (e.g. aa, g) is done via the sanitize_species
        method from the util.species module.
        """
        # tokenize the rate string and apply any sanitizations
        tokenized = map(sanitize_species,
                        rate_breaker.match(self.rxnString).groups())
        reactants = tokenized[:2]
        products = tokenized[2:]
        print 'reactants', reactants
        print 'products', products
        # target = sanitize_species(self.rxnString.split('(')[0])
        # internals = self.rxnString.split('(')[1].split(')')[0].split(',')
        # reactants, products = [sanitize_species(spec) for spec in internals]
        # endState = sanitize_species(self.rxnString.split(')')[1])
        

        self.reactants = reactants  #.split(PLUS)
        self.products = products  #.split(PLUS)

        # # this is of the form a + b -> c + d, with special character handling
        # rateString = form_rate_string(*tokenized)

        # # split this into reactants and products, for easier searching of the
        # # various permutations in the ReacLib database
        # self.reactants = rateString.split(rl.TO)[0].split(rl.PLUS)
        # self.products = rateString.split(rl.TO)[1].split(rl.PLUS)


    def _build_rxn_data(self):
        # find this reaction in the reaction rate file
        rate_data = webnucleo.get_rate_data(self)
        # rate data is stored in several formats
        rate_builders = {"non_smoker_fit": webnucleo.build_non_smoker_rate,
                         "single_rate": webnucleo.build_single_rate,
                         "rate_table": webnucleo.build_rate_table_rate}
        for rate_type, rate_builder in rate_builders.iteritems():
            this_rate = rate_data.find(rate_type)
            # if the storage type matches, then read the data and build
            # the rate function, a function of temperature only
            if this_rate:
                rate_builder(self, rate_data)

        # we need to build the reverse factor; for a non-reverse rate, this is
        # just unity
        # this is a function of density and temperature
        self.reverse_factor = webnucleo.build_reverse_rate_function(self)

        # the full reaction rate
        def _full_rate(reaction, temperature, density):
            return (reaction.reverse_factor(temperature, density) *
                    reaction.forward_rate(temperature))
        self.rate = _full_rate

    def _build_qvalue(self):
        """
        Calculate the Q-value of the reaction, including modifiers for 
        electron-capture and beta+ decay.
        """
        qvalue = np.sum([isotope.mass_excess
                         for isotope in self.isotope_reactants])
        qvaue -= np.sum([isotope.mass_excess
                         for isotope in self.isotope_products])
        # if this is a beta decay, we lose twice electron mass
        if self.is_betaplus:
            qvalue -= 2 * electron_mass * light_speed * light_speed
        # if this is electron capture, we gain twice electron mass
        elif self.is_electron_capture:
            qvalue += 2 * electron_mass * light_speed * light_speed
        self.qvalue = qvalue

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


# class ReacLibReaction(Reaction):
#     """
#     For ReacLib[1] reaction rate objects.  These will be stored in an XML
#     file that we only want to parse once, so we pass just the name of the
#     isotope -- NOT the Isotope object -- to the baseclass.

#     [1] https://groups.nscl.msu.edu/jina/reaclib/db/
#     """
#     def __init__(self, nucRxnString):
#         """
#         nucRxnString is a reaction descriptor in the form of a(b,c)d
#         """
#         self.nucRxnString = nucRxnString
#         # build target, reactants and products
#         self._parse_rate_string()

#         # just the nuclides, ignoring leptons and gammas
#         nuclide_reactants = [reactant for reactant in self.reactants
#                              if reactant not in ['gamma', 'lepton']]
#         nuclide_products = [product for product in self.products
#                             if product not in ['gamma', 'lepton']]

#         # separate into unique values
#         isotopes_cnt = collections.Counter(nuclide_reactants +
#                                            nuclide_products)
#         uniq_reactants = set(nuclide_reactants)
#         uniq_products = set(nuclide_products)
#         species = [(-1*isotopes_cnt[isotope], isotope)
#                    for isotope in uniq_reactants]
#         species += [(isotopes_cnt[isotope], isotope)
#                     for isotope in uniq_products]

#         super(ReacLibReaction, self).__init__(species)


#     def build_rxn_rate(self, xml_root):
#         """
#         xml_root is the root node of the XML document containing all the
#         reaction rate information in the ReacLib database.
#         """
#         # xpath_str = '//reaction[%s]'
#         # reactants = ['reactant="%s"' % reactant.lower()
#         #              for reactant in self.reactants]
#         # products = ['product="%s"' % product.lower()
#         #             for product in self.products]
#         # xpath_str = xpath_str % ' and '.join(reactants + products)
#         # search for reactions that have these specific reactants and products
#         possible_reactions = xml_root.xpath(_build_xpath(self.reactants,
#                                                          self.products))
#         # if possible_reactions is an empty list, then this might be
#         # a reverse rate, which isn't stored in the database.  need to
#         # construct the reverse rate, if possible
#         if not possible_reactions:
#             possible_reactions = xml_root.path(_build_xpath(self.products,
#                                                             self.reactants))
#             if not possible_reactions:
#                 errStr = "%s\n or it's reverse not found in ReacLib database."
#                 raise RuntimeError(errStr % str(self))

#             # FIXME - implement this
#             _build_reverse()

#         # possible_reactions is NOT a list of reactions with ONLY these
#         # products and reactants, so we need to filter some out by making
#         # sure the NUMBER of products and reactants match, respectively
#         rxns_same_no_reactants = [rxn for rxn in possible_reactions if
#                                   int(rxn.xpath('count(.//reactant)')) ==
#                                   len(self.reactants)]
#         this_reaction = [rxn for rxn in rxns_same_no_reactants if
#                          int(rxn.xpath('count(.//product)')) ==
#                          len(self.products)]
#         if len(this_reaction) != 1:
#             errStr = "Found multiple reactions for \n%s\n" % str(self)
#             errStr += xpath_str
#             raise RuntimeError(errStr)
#         this_reaction = this_reaction[0]
#         # now read in the recommended fit parameters
#         aFactors = {}
#         for a_param in ['a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7']:
#             aFactors[a_param] = [float(a_val.text) for a_val in
#                                  this_reaction.xpath('//%s' % a_param)]

#         self.rate = rl.build_rate_function(aFactors)

#     def reverse_factor(self, t, rho):
#         return 1.0


# def _build_xpath(reactants, products):
#     # use XPath syntax to search for this specific rate based on
#     # the reactants and products
#     xpath_str = '//reaction[%s]'
#     reacts = ['reactant="%s"' % reactant.lower()
#               for reactant in reactants]
#     prods = ['product="%s"' % product.lower()
#              for product in products]
#     return xpath_str % ' and '.join(reacts + prods)

# def _build_reverse():
#     """
#     Need to implement reverse rates here
#     """
#     pass
