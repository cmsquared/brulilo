"""
This module provides an interface for acquiring reaction rates from the
ReacLib[1] database.  It queries the database over HTTP, in a somewhat
inefficient fashion.

[1] https://groups.nscl.msu.edu/jina/reaclib/db/
"""
import urllib
import sys
import re
import itertools

# ReacLib format 1 - this is all we support currently
RL_format = 'cfreaclib'

# some oft-used constants
PLUS = " + "
TO = " -> "
MAX_PER_PAGE = 9999

# these are used to parse HTML files and species names
specZAFinder = re.compile(r'(\D+)(\d+)')
rateFinderString = '<td><a href="(.*)">%s'
idAndFNFinder = re.compile(r"'&rateID=(\d+)&filename=(.*)'"
                           r"\+ \(this.action.value=='xml'")

# periodic table, somewhat broken up appropriately
isotope_lut = ['H', 'He',
               'Li', 'Be',
               'B', 'C', 'N', 'O', 'F', 'Ne',
               'Na', 'Mg',
               'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
               'K', 'Ca',
               'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co',  # break
               'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
               'Rb', 'Sr',
               'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh',  # break
               'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
               'Cs', 'Ba',
               # Lanthanide Series
               'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu',  # break
               'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
               'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt',  # break
               'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
               'Fr', 'Ra',
               # Actinide Series
               'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am',  # break
               'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr',
               'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',  # break
               'Rg', 'Cn', 'Uut', 'Fl', 'Uup', 'Lv', 'Uus', 'Uuo']
# this is useful for going from a species name to a Z value
Zdict = {}
for i, species in enumerate(isotope_lut):
    Zdict[species] = i+1

# ReacLib uses some common symbols for proton, deuteron, etc.
specialCharZA = {'p': (1, 1),
                 'd': (1, 2),
                 't': (1, 3),
                 'a': (2, 4),
                 'n': (0, 1)}


def get_Z_A(spec):
    """
    Parses a species like 'He4' or 'a' to get the A and Z for searching.
    """
    foundList = re.findall(specZAFinder, spec)
    # if it is empty, then we have some special chars
    if not foundList:
        return specialCharZA[spec]
    else:
        spec, A = foundList[0]
        Z = Zdict[spec]
    return Z, int(A)


def sanitize_species(speciesString):
    """
    Takes a 'species' from a reaction string and parses it properly.
    e.g:
      'nn' really means 'n + n', so that is what will be returned.
      'g'  really means a photon, which is omitted
      'He4' is a bonafide isotope, so just return it
    """
    # ignore photons and empty strings
    if speciesString.strip() == '' or speciesString == 'g':
        return ''
    # if there are capitals in the string, then assume this is a bonafide
    # isotope
    if speciesString.lower() != speciesString:
        return speciesString
    # now, if there are any 'g's left, set them to empty
    speciesString = speciesString.replace('g', '')
    # if it is a single species left, correct alphas and return
    if len(speciesString) is 1:
        return speciesString.replace('a', 'He4')
    else:
        # multiple species
        ret = []
        for spec in speciesString:
            ret.append(spec.replace('a', 'He4'))
        ret = PLUS.join(ret)
        return ret


def form_rate_string(target, reactants='', products='', endState=''):
    """
    Takes the various parts of a reaction and forms the stochiometric version,
    i.e.
      target + reactants -> endState + products
    """
    rstring = [target, ]
    if reactants is not '':
        rstring += [PLUS, reactants]
    rstring.append(TO)
    if endState is not '':
        rstring.append(endState)
    if products is not '':
        rstring += [PLUS, products]

    # this is of the form a + b -> c + d, with special character handling
    rstring = ''.join(rstring)
    return rstring


def search_for_reactions_involving_target(Z, A):
    """
    This will call the ReacLib search functionality specifying the target
    isotope as the (Z,A) isotope.  This should return a page with all such
    reactions.
    """
    baseURL = 'https://groups.nscl.msu.edu/jina/reaclib/db/results.php'
    data = {'lowz[]': str(Z), 'highz[]': str(Z),
            'lowm[]': str(A), 'highm[]': str(A),
            'perPage': MAX_PER_PAGE}
    data = urllib.urlencode(data)
    try:
        req = urllib.urlopen(baseURL, data)
    except IOError:
        print ("Couldn't open %s with the following data \n %s" %
               (baseURL, data))
        sys.exit()

    return req.read()
