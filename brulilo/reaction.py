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
    For ReacLib[1] reaction rate objects. This stores, parses and searches for
    rate information over HTTP from the ReacLib database.

    [1] https://groups.nscl.msu.edu/jina/reaclib/db/
    """
    def __init__(self, nucRxnString):
        """
        nucRxnString is a reaction descriptor in the form of a(b,c)d
        """
        self.nucRxnString = nucRxnString
        # build target, reactants and products
        self._parse_rate_string()

        # here we need to build the intervening Isotopes to pass to the base
        # init method; we assume no Isotopes are catalysts
        all_reactants = [Isotope(spec) for spec in self.reactants]
        all_products = [Isotope(spec) for spec in self.products]
        isotopes_cnt = collections.Counter(map(str,
                                               all_reactants + all_products))
        uniq_reactants = set(all_reactants)
        uniq_products = set(all_products)
        uniq_isotopes = uniq_reactants | uniq_products
        print self.nucRxnString
        print 'isotopes_cnt', isotopes_cnt
        print 'uniq_isotopes', map(str, uniq_isotopes)
        # isotopes = [Isotope(spec) for spec in uniq_isotopes]
        # # assuming there aren't any catalytic species in the reaction
        # reactants = [isotope for isotope in isotopes
        #              if str(isotope) in self.reactants]
        # products = [isotope for isotope in isotopes
        #             if str(isotope) in self.products]
        species = [(-1*isotopes_cnt[str(isotope)], isotope)
                   for isotope in uniq_reactants]
        species += [(isotopes_cnt[str(isotope)], isotope)
                    for isotope in uniq_products]
        print 'species'
        for spec in species:
            print spec[0], str(spec[1])
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
        target = self.nucRxnString.split('(')[0]
        internals = self.nucRxnString.split('(')[1].split(')')[0].split(',')
        reactants, products = [rl.sanitize_species(spec) for spec in internals]
        endState = self.nucRxnString.split(')')[1]

        # this is of the form a + b -> c + d, with special character handling
        rateString = rl.form_rate_string(target, reactants, products, endState)

        # split this into reactants and products, for easier searching of the
        # various permutations in the ReacLib database
        self.reactants = rateString.split(rl.TO)[0].split(rl.PLUS)
        self.products = rateString.split(rl.TO)[1].split(rl.PLUS)

        self.target = target

    def _find_specific_reaction(self,pageContainingAllRatesForTarget):
        """
        There is not a standard way reactions are written in the pages returned
        from a search query; i.e. n + n + He4 -> He6 could also be written as
        n + He4 + n -> He6 or He4 + n + n -> He6.
        
        This method searches for all the possible rates given some target
        nucleus, and finds the correct rate by looking through all possible
        permutations of reactants and products.  If there are multiple matches
        (sometimes the case with weak rates), then the user is queried about
        which one they prefer.
        """
        reactantsPerms = set(itertools.permutations(self.reactants))
        productsPerms = set(itertools.permutations(self.products))
        allPerms = itertools.product(reactantsPerms,productsPerms)

        # loop over all the permutations and find a matching string
        # we store the reaction rate URLs in totalMatches
        totalMatches = []
        for rperm,pperm in allPerms:
            target = rperm[0]
            reactants  = ''
            if len(rperm) > 1: reactants = PLUS.join(rperm[1:])
            products = ''
            if len(pperm) > 1: products = PLUS.join(pperm[:-2])
            endState = pperm[-1]
            rstring = form_rate_string(target,reactants,products,endState)

            matches = re.findall(rateFinderString % rstring.replace('+','\+'),
                                 pageContainingAllRatesForTarget)
            if matches: totalMatches += matches

        if not totalMatches: 
            print "Error, couldn't find reaction for %s" % self.nucRxnString
            sys.exit()

        thisMatch = totalMatches[0]

        # if there are multiple choices, then prompt the user for a preference
        # this happens with weak decays that can be either beta+ or e capture
        if len(totalMatches) > 1:
            print "\nThere are multiple matches for the rxn you specified."
            print "Note that this is likely the case for weak rates that can"
            print "either beta^+ decay or electron capture."
            print "\nWhich rate would you prefer?\n"
            for i,m in enumerate(totalMatches):
                print '%s) %s' % (i+1,m)
            iselect = raw_input()
            thisMatch = totalMatches[i-1]
        
        self.rateURL = thisMatch

    def _get_ID_filename(self):
        """
        Given a reaction rate URL, this function will parse the site for the 
        reaction rate ID and the filename in the database.  This can then be
        used to search for the exact rate file.
        """
        data = urllib.urlopen(self.rateURL).read()
        matches = idAndFNFinder.findall(data)
        if not matches:
            print "Couldn't find ID and filename from %s" % URL
        self.rateID, self.filename = matches[0]

    def get_RL_file(self,pageContainingAllRatesForTarget):
        """
        Given an HTML page containing all the reactions associated with this
        particular Reaction's target, find this particular reaction and 
        download the ReacLib rate file containing the fitting parameters.
        """
        # get the rate URL for this particular Reaction
        self._find_specific_reaction(pageContainingAllRatesForTarget)
        # parse the rate URL for the ReacLib ID and filename
        self._get_ID_filename()

        # download and write it to a file
        baseURL = ('https://groups.nscl.msu.edu/jina/reaclib/db/' +
                   'difout.php?action=%s' +
                   '&rateID=%s&filename=%s&no910=0')
        fullURL = baseURL % (RL_format, self.rateID, self.filename)
        req = urllib.urlopen(fullURL)
        try:
            fh = open(self.filename,'w')
        except IOError:
            print "Couldn't open %s for writing!" % self.filename
            sys.exit()
        fh.write(req.read())
        fh.close()
