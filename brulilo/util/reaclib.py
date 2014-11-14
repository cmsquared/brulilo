"""
This module provides an interface for acquiring reaction rates from the
ReacLib[1] database.  Currently, we only consider data in the XML format 
described by Webnucleo[2].

[1] https://groups.nscl.msu.edu/jina/reaclib/db/
[2] http://nucleo.ces.clemson.edu/
"""
import re
import os.path
import lxml.etree as etree
import numpy as np

import brulilo
from constants import MeV2erg
from progressbar import IntProgressBar








def get_nuclear_data(network):
    """
    Grab all the nuclear data for all Isotopes in the Network.
    """
    fn = os.path.join(os.path.dirname(brulilo.__file__),
                      nuc_data_file)
    nuc_data_root = etree.parse(fn)
    pbar = IntProgressBar("Getting nuclear data", len(network.isotopes))
    for isotope in network.isotopes:
        # find this element in the XML file based on its Z and A
        xpath_str = "//nuclide[z=%d and a=%d]" % (isotope.Z, isotope.A)
        this_element = nuc_data_root.xpath(xpath_str)
        if len(this_element) != 1:
            errStr = "Didn't find a single isotope for (%d,%d)\n" % (isotope.Z,
                                                                     isotope.A)
            errStr += "Found" + ' '.join(this_element)
            raise RuntimeError(errStr)
        this_element = this_element[0]
        # get the mass excess
        isotope.mass_excess = (float(this_element.find("mass_excess").text) *
                               MeV2erg)
        # we need the following for calculating reverse rates
        # get the spin
        isotope.spin = float(this_element.find("spin").text)
        # get the partition table for interpolation
        # this is going to be an array of t9, log10(f) entries
        # -- see Webnucleo documentation for more info
        partition_table = this_element.find("partf_table")
        part_t9 = [float(t9.text) for t9 in partition_table.xpath("point/t9")]
        part_lf = [float(lf.text) for lf in
                   partition_table.xpath("point/log10_partf")]
        isotope.partition_table = np.array(zip(part_t9, part_lf),
                                           dtype='float64')
        pbar.update(str(isotope))

def get_rate_data(network):
    """
    Grab all the nuclear reaction data for all Reactions in the Network.
    Here we also build the rate functions for each Reaction, including
    the machinery for reverse rates.
    """
    fn = os.path.join(os.path.dirname(brulilo.__file__),
                      rxn_data_file)
    rxn_data_root = etree.parse(fn)
    for reaction in network.reactions:
        print reaction
        # xpath = "//reaction[%s]"
        # # first add the reactant nuclides + gammas
        # xpath_specific = AND.join(["reactant='%s'" % reactant
        #                            for reactant in reaction.reactants
        #                            if reactant != "lepton"])
        # # now the product nuclides + gammas
        # xpath_specific += AND + AND.join(["product='%s'" % product
        #                                   for product in reaction.products
        #                                   if product != "lepton"])
        # # 'lepton' indicates that it is something from leptons, so search all
        # # reactants first
        # if 'lepton' in reaction.reactants:
        #     xpath_specific += AND + "(" + OR.join(["reactant='%s'" % lepton
        #                                            for lepton in leptons]) + \
        #                                                ")"
        # if 'lepton' in reaction.products:
        #     xpath_specific += AND + "(" + OR.join(["product='%s'" % lepton
        #                                            for lepton in leptons]) + \
        #                                                ")"
        search_string = _build_rate_xpath_str(reaction.reactants,
                                              reaction.products)
        possible_reactions = rxn_data_root.xpath(search_string)
        npossible = len(possible_reactions)
        this_reaction = None
        print search_string
        # if we didn't find anything, this is possibly a reverse rate
        reaction.is_reverse = False
        if npossible == 0:
            reaction.is_reverse = True
            this_reaction = _try_to_grab_reverse(reaction, rxn_data_root)
        elif npossible == 1:
            this_reaction = possible_reactions[0]
        # if we found two rates, then this is probably an ambiguous
        # case that needs to check the qualifier
        # TODO -- cleanup
        elif npossible == 2:
            this_reaction = _address_qualifier(possible_reactions)
        else:
            errStr = ("Found more than two reactions for: %s\n" %
                      reaction.nucRxnString)
            for rxn in possible_reactions:
                errStr += "\n"
                for child in rxn:
                    errStr += "%s %s" % (str(child), child.text)
            raise RuntimeError(errStr)
        _build_rate_function(reaction, this_reaction)


def _build_rate_function(rxn, xml_rxn):
    """
    Here we parse the XML reaction object to find either a single_rate tag,
    a rate_table tag, or a non_smoker_fit.  We build the appropriate rate
    function and add it to the Reaction object along with the reaction's 
    Q-value.  Special handling is done for reverse rates, and Q-values are
    modified based on if we are a beta+ rate or electron capture.
    """
    rate_builders = {"non_smoker_fit": _build_non_smoker_rate,
                     "single_rate": _build_single_rate,
                     "rate_table": _build_rate_table_rate}
    print 'xml_rxn', xml_rxn
    for rate_type, builder in rate_builders.iteritems():
        this_rate = xml_rxn.xpath(rate_type)
        if this_rate:
            rxn.rate = builder(rxn, this_rate[0])
            return

def _build_non_smoker_rate(rxn, xml_rxn):
    """
    Read in the 'a' factors to the non-smoker fit and build the rate function.
    TODO: This currently doesn't honor the lower and upper temperature bounds..
    """
    aList = ['a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7']
    aFactors = {}
    for a_param in aList:
        aFactors[a_param] = [float(a_val.text) for a_val in
                             xml_rxn.xpath('.//%s' % a_param)]
    aFacs = np.array(zip(*[aFactors[a] for a in aList]))
    def _rate_function(self, temperature, density):
        t9 = temperature / 1e9
        t9i = 1./t9
        tfactors = np.array([1.0,
                             t9i, t9i**(1./3.),
                             t9**(1./3.), t9, t9**(5./3.),
                             np.log10(t9)], dtype='float64')
        rate = np.exp(aFacs * tfactors)
        rate = np.sum(rate)
        if self.is_reverse:
            rate *= self.reverse_factor(temperature, density)
        return rate
    return _rate_function

def _build_single_rate(rxn, xml_rxn):
    """
    Read in the single rate and build the function.  A function is overkill
    in this case, but keeps the API consistent with the other formats.
    """
    single_rate = xml_rxn.xpath(".//single_rate")
    if not single_rate:
        errString = "Couldn't find the single_rate for %s" % str(rxn)
        raise RuntimeError(errString)
    single_rate = float(single_rate[0])
    def _rate_function(self, temperature, density):
        if self.is_reverse:
            single_rate *= self.reverse_factor(temperature, density)
        return single_rate
    return _rate_function

def _build_rate_table_rate(rxn, xml_rxn):
    raise NotImplementedError


def _build_rate_xpath_str(lhs, rhs):
    """
    this builds an XPath string for searching through the XML using the
    lhs as the reactants and the rhs as the products
    """
    xpath = "//reaction[%s]"
    # first add the reactant nuclides + gammas
    xpath_specific = AND.join(["reactant='%s'" % reactant
                               for reactant in lhs if reactant != "lepton"])
    # now the product nuclides + gammas
    xpath_specific += AND + AND.join(["product='%s'" % product
                                      for product in rhs
                                      if product != "lepton"])
    # 'lepton' indicates that it is something from leptons, so search all
    # reactants first
    if 'lepton' in rhs:
        xpath_specific += AND + "(" + OR.join(["reactant='%s'" % lepton
                                               for lepton in leptons]) + \
                                                   ")"
    if 'lepton' in lhs:
        xpath_specific += AND + "(" + OR.join(["product='%s'" % lepton
                                               for lepton in leptons]) + \
                                                   ")"

    return xpath % xpath_specific.lower()
    

def _try_to_grab_reverse(rxn, xml_root):
    # reverse the reactants and products
    reverse_search_str = _build_rate_xpath_str(rxn.products, rxn.reactants)
    reverse_reactions = xml_root.xpath(reverse_search_str)
    if not reverse_reactions:
        errString = "Could not find forward or reverse rate for: %s" % str(rxn)
        raise RuntimeError(errString)
    if len(reverse_reactions) > 1:
        reverse_reaction = _address_qualifier(reverse_reactions)
    else:
        reverse_reaction = reverse_reactions[0]
    return reverse_reaction

def _address_qualifier(reactions):
    pass
