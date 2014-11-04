"""
This module provides an interface for acquiring reaction rates from the
ReacLib[1] database.

[1] https://groups.nscl.msu.edu/jina/reaclib/db/
"""
import re

# these files store data about the rates and nuclei in the ReacLib database
nuc_data_file = "data/webnucleo_nuc_v2.0.xml"
rxn_data_file = "data/20141031default.webnucleo.xml"

# some oft-used constants
PLUS = " + "
TO = " -> "

# this is used to parse species names
specZAFinder = re.compile(r'(\D+)(\d+)')

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
    if speciesString.strip() in ['', 'g']:
        return ''
    # if there are capitals in the string, then assume this is a bonafide
    # isotope
    if speciesString.lower() != speciesString:
        return speciesString
    # now, if there are any 'g's left, set them to empty
    speciesString = speciesString.replace('g', '')
    # if it is a single species left, correct specials and return
    if len(speciesString) == 1:
        return _fix_special_species(speciesString)

    # multiple species; correct specials
    ret = []
    for spec in speciesString:
        ret.append(_fix_special_species(spec))
    ret = PLUS.join(ret)
    return ret

def _fix_special_species(species_string):
    """
    Fix the cases of special characters; i.e. t --> H3, a --> He4, etc.
    """
    # don't do anything to n
    if species_string == 'n':
        return species_string

    if species_string in specialCharZA:
        Z, A = specialCharZA[species_string]
        return "%s%d" % (isotope_lut[Z-1], A)

    return species_string

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

def build_rate_function(aFactors):
    """
    aFactors is a dictionary whose keys are the name of the fit parameters
    (e.g. 'a1' or 'a5') and the values are lists of the fit parameter values,
    one entry for each set in the ReacLib rate
    """
    aFacs = np.array(zip(aFactors['a1'], aFactors['a2'], aFactors['a3'],
                         aFactors['a4'], aFactors['a5'], aFactors['a6'],
                         aFactors['a7']))
    def _rate_function(temperature):
        t9 = temperature / 1e9
        t9i = 1./t9
        tfactors = np.array([1.0,
                             t9i, t9i**(1./3.),
                             t9**(1./3.), t9, t9**(5./3.),
                             np.log10(t9)], dtype='float64')
        rate = np.exp(aFacs * tfactors)
        return np.sum(rate)
        
    return _rate_function
