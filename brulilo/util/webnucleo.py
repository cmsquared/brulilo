"""
This provides some parsing and handling of the Webnucleo[1] data format.

[1] http://nucleo.ces.clemson.edu/
"""
import lxml.etree as etree
from scipy.interpolate import interp1d
import numpy as np
import os.path

import brulilo
from .constants import avogadro, light_speed, boltzmann, planck_bar, amu


class WebnucleoDataParser(object):
    __base_dir = os.path.dirname(brulilo.__file__)

    _nuc_data_file = os.path.join(__base_dir,
                                  "data/webnucleo_nuc_v2.0.xml")
    _rxn_data_file = os.path.join(__base_dir,
                                  "data/20141031default.webnucleo.xml")

    # cache the parsing of the data files
    _nuc_xml_root = None

    @property
    def nuc_data_file(self):
        if self._nuc_xml_root is None:
            self._nuc_xml_root = etree.parse(self._nuc_data_file)
        return self._nuc_xml_root

    _rxn_xml_root = None

    @property
    def rxn_data_file(self):
        if self._rxn_xml_root is None:
            self._rxn_xml_root = etree.parse(self._rxn_data_file)
        return self._nuc_xml_root

    # these are used in calculating binding energies
    _proton_mass_excess = None

    @property
    def proton_mass_excess(self):
        if self._proton_mass_excess is None:
            xpath_str = "//nuclide[z=1 and a=1]"
            proton = self.nuc_data_file.xpath(xpath_str)[0]
            self._proton_mass_excess = float(proton.find("mass_excess").text)
        return self._proton_mass_excess

    _neutron_mass_excess = None

    @property
    def neutron_mass_excess(self):
        if self._neutron_mass_excess is None:
            xpath_str = "//nuclide[z=0 and a=1]"
            neutron = self.nuc_data_file.xpath(xpath_str)[0]
            self._neutron_mass_excess = float(neutron.find("mass_excess").text)
        return self._neutron_mass_excess

    def get_isotope_data(self, isotope):
        """
        Finds and returns a specific isotope within the nuclide data file.
        The returned object is an etree.Element object.
        """
        # find this element in the XML file
        xpath_str = "//nuclide[z=%d and a=%d]" % (isotope.Z, isotope.A)
        this_isotope = self.nuc_data_file.xpath(xpath_str)
        # if we didn't find a single entry, then this is an error
        if len(this_isotope) != 1:
            errString = ("Didn't find a proper entry for isotope %s" %
                         str(isotope))
            raise RuntimeError(errString)
        return this_isotope[0]

    def get_rate_data(self, reaction):
        """
        Finds and returns a specific reaction rate within the data file.
        The returned object is an etree.Element object.
        """
        # find this reaction in the XML file
        xpath_str = _build_rxn_xpath_search(reaction.reactants,
                                            reaction.products)
        this_reaction = self.rxn_data_file.xpath(xpath_str)
        # if we didn't find anything, then this is a reverse rate
        if not this_reaction:
            reaction.is_reverse = True
            # swap the reactants and products and re-search
            xpath_str = _build_rxn_xpath_search(reaction.products,
                                                reaction.reactants)
            this_reaction = self.rxn_data_file.xpath(xpath_str)
            # now if THIS is empty, we have an error
            if not this_reaction:
                errString = "Couldn't find either a forward or reverse"
                errString += " rate for\n %s" % reaction.rxnString
                raise RuntimeError(errString)
        return this_reaction[0]

    def build_non_smoker_rate(self, reaction, reaction_xml):
        """
        This reaction's data is stored as a non-smoker fit in the 
        reaction_xml data.  The units are rate per interaction pair or 
        multiplet per second.  Parse it and build the forward rate function.
        """
        # TODO -- honor min/max temperatures
        aList = ["a%s" % i+1 for i in range(7)]
        aFactors = {}
        for a in aList:
            aFactors[a] = map(float, [a_val.text for a_val in
                                      reaction_xml.xpath('.//%s' % a)])
        aFacs = np.array(zip(*[aFactors[a] for a in aList]))
        
        def _forward_rate_function(reaction, temperature):
            t9 = temperature / 1e9
            t9i = 1./t9
            tfactors = np.array([1.0,
                                 t9i, t9i**(1./3.),
                                 t9**(1./3.), t9, t9**(5./3.),
                                 np.log10(t9)], dtype='float64')
            rate = np.exp(aFacs * tfactors)
            return np.sum(rate)
        reaction.forward_rate = _forward_rate_function

    def build_single_rate(self, reaction, reaction_xml):
        """
        This reaction's data is stored as a single rate in the reaction_xml
        data.  The units are rate per nuclide per second. Parse it and build 
        the forward rate function.
        """
        single_rate = float(reaction_xml.xpath(".//single_rate")[0].text)
        def _forward_rate_function(reaction, temperature):
            return single_rate
        reaction.forward_rate = _forward_rate_function

    def build_rate_table_rate(self, reaction, reaction_xml):
        """
        This reaction's data is stored in a table of (t9, rate) pairs.  The
        units are rate per interaction pair or multiplet per second.  Parse
        it and build the forward rate function.
        """
        rtable = reaction_xml.find("rate_table")
        rt9 = map(float, [t9.text for t9 in rtable.xpath("point/t9")])
        rrate = map(np.log10, [rate.text
                               for rate in rtable.xpath("point/rate")])
        # stellar enhancement factor; accounts for excited states
        rsef = map(np.log10, [sef.text
                              for sef in rtable.xpath("point/sef")])
        # this is the total rate
        rtotal = [rate + sef for rate, sef in zip(rrate, rsef)]
        mint9 = min(rt9)
        maxt9 = max(rt9)
        # now the fit
        rfit = interp1d(rt9, rtotal, kind='cubic')
        def _forward_rate_function(reaction, temperature):
            t9 = temperature / 1e9
            # no extrapolation
            t9 = min(maxt9, max(mint9, t9))
            return 10**rfit(t9)
        reaction.forward_rate = _forward_rate_function

    def build_reverse_rate_function(self, reaction):
        """
        Build a function that calculates the reverse rate factor from detailed
        balance.  This factor is multiplied by the forward rate to get the
        total rate.
        """
        # first, if this isn't a reverse reaction, then the multiplicative
        # factor is just unity
        if not reaction.is_reverse:
            return lambda temperature, density: 1.0
        # if the forward reaction is weak, then this is not reversible
        if reaction.is_weak:
            return lambda temperature, density: 0.0
        # bonafide, non-weak reverse reaction
        def _reverse_factor(rxn, temperature, density):
            # common factors
            factor1 = 1. / (boltzmann * temperature)
            factor2 = (1. /
                       (factor1 * 2 *  np.pi *
                        (planck_bar * light_speed)**2))**(3./2.)
            factor3 = 1. / (avogadro * density)
            # density weighting
            dexp = (len(rxn.isotope_reactants)
                    - len(rxn.isotope_products)) * np.log(density)
            # reactants
            dexp += np.sum([np.log(iso.partition_function(temperature) *
                                   factor2 * factor3 *
                                   (iso.A * amu * light_speed**2 +
                                    iso.mass_excess)**(3./2.))
                            + iso.binding_energy * factor1 for iso in
                            rxn.isotope_reactants])
            # products
            dexp -= np.sum([np.log(iso.partition_function(temperature) *
                                   factor2 * factor3 *
                                   (iso.A * amu * light_speed**2 +
                                    iso.mass_excess)**(3./2.))
                            + iso.binding_energy * factor1 for iso in
                            rxn.isotope_products])
            # TODO -- account for duplicate reactants/products
            return np.exp(dexp)

        return _reverse_factor

# instantiate a class to handle the parsing
webnucleo = WebnucleoDataParser()


def _build_rxn_xpath_search(reactants, products):
    """
    Helper function for searching for forward/reverse rates within an XML file.
    """
    xpath_str = "//reaction[%s]"
    reacts = ["reactant='%s'" % reactant.lower() for reactant in reactants]
    prods = ["product='%s'" % product.lower() for product in products]
    return xpath_str % ' and '.join(reacts + prods)
