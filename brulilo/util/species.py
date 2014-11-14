"""
Some utilities and maps for handling the various species we encounter - both
isotopes and leptons and photons.
"""
import re


# this dictionary is a map between the syntax in our reaction rate file
# and what Webnucleo uses
rxn_to_WN_map = {'e-': 'electron', 'e+': 'positron',
                 'nu_e': 'neutrino_e', 'nu_e_bar': 'anti-neutrino_e',
                 'mu-': 'mu', 'mu+': 'anti-mu',
                 'nu_mu': 'neutrino_mu', 'nu_mu_bar': 'anti-neutrino_mu',
                 'tau-': 'tau', 'tau+': 'anti-tau',
                 'nu_tau': 'neutrino_tau', 'nu_tau_bar': 'anti-neutrino_tau',
                 'g': 'gamma'}
non_nuclides = rxn_to_WN_map.values()
leptons = [l for l in non_nuclides if l != "gamma"]

# periodic table, somewhat broken up appropriately
element_lut = ['n',
               'H', 'He',
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
# these correspond to the minimum and maximum A for each element in element_lut
# some of these aren't really stable enough for reactions, but oh well...
# labels are the last element in that line
# data taken from webnucleo file
isotope_A_ranges = [
    (1, 1),  # n
    (1, 3), (3, 10),  # He
    (6, 11), (7, 14),  # Be
    (8, 19), (9, 22), (11, 23), (13, 30), (14, 37), (16, 41),  # Ne
    (18, 44), (19, 47),  # Mg
    (21, 51), (22, 54), (23, 57), (24, 60), (25, 63), (27, 67),  # Ar
    (29, 70), (30, 73),  # Ca
    (32, 76), (34, 80), (36, 83), (38, 86), (40, 89), (42, 92), (44, 96),  # Co
    (46, 99), (48, 102), (51, 105), (53, 108), (55, 112), (57, 115),  # As
    (59, 118), (61, 121), (63, 124),  # Kr
    (66, 128), (68, 131),  # Sr
    (70, 134), (72, 137), (74, 140), (77, 144), (79, 147), (81, 150),  # Ru
    (83, 153),  # Rh
    (86, 156), (88, 160), (90, 163), (92, 166), (94, 169), (97, 172),  # Sb
    (99, 176), (101, 179), (103, 182),  # Xe
    (106, 185), (108, 189),  # Ba
    (110, 192), (113, 195), (115, 198), (118, 201), (120, 205),  # Pm
    (123, 208), (125, 211), (128, 214), (130, 218), (133, 221),  # Dy
    (136, 224), (138, 227), (141, 230), (143, 234), (146, 237),  # Lu
    (149, 240), (151, 243), (154, 247), (156, 250), (159, 253),  # Os
    (162, 256), (165, 260), (167, 263), (170, 266), (173, 269),  # Tl
    (175, 273), (178, 276), (181, 276), (184, 279), (188, 269),  # Rn
    (191, 280), (194, 283),  # Ra
    (197, 288), (200, 293), (213, 296), (205, 295), (221, 302),  # Np
    (211, 303), (216, 309), (217, 312), (226, 315), (223, 319),  # Cf
    (226, 322), (228, 325), (240, 328), (234, 332), (241, 335),  # Lr
    (240, 337), (243, 337), (244, 337), (247, 337), (250, 337),  # Hs
    (253, 337), (256, 337), (259, 337), (262, 337), (268, 337),  # Uut
    (270, 337), (273, 337), (276, 337), (280, 337), (283, 337)   # Uuo
    ]
# a list of all the possible isotopes
isotope_lut = [element_lut[0], ]
for i, A_range in enumerate(isotope_A_ranges):
    if i == 0:  # skip neutron
        continue
    isotope_lut.extend("%s%d" % (element_lut[i], A)
                       for A in range(A_range[0], A_range[1]+1))

# this is useful for going from a species name to a Z value
Zdict = {}
for i, species in enumerate(element_lut):
    Zdict[species] = i

# we use some common symbols for proton, deuteron, etc.
specialCharZA = {'p': (1, 1),
                 'd': (1, 2),
                 't': (1, 3),
                 'a': (2, 4),
                 'n': (0, 1)}

# these are used to parse species names
_specZAFinder = re.compile(r'(\D+)(\d+)')
_specSplitter = re.compile(r'([A-Z][^A-Z]*)')

# some oft-used constants
PLUS = " + "
TO = " -> "
AND = " and "
OR = " or "


def get_Z_A(spec):
    """
    Parses a species like 'He4' or 'a' to get the A and Z for searching.
    """
    foundList = re.findall(_specZAFinder, spec)
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
      'nn' really means 'n + n', so that is what will be returned
      'g'  really means a photon, so return 'gamma'
      'He4' is a bonafide isotope, so just return it
      'e+' really means there was a positron
      'aO20' really means 'He4 + O20'
    """
    # TODO -- optimize/clarify this

    # check for pure species
    if speciesString in isotope_lut:
        return speciesString
    # check for pure non-isotopes, converting to Webnucleo syntax
    if speciesString in rxn_to_WN_map:
        return rxn_to_WN_map[speciesString]
    # check for pure special characters; this doubly checks neutron...
    if speciesString in specialCharZA:
        return _fix_special_species(speciesString)
    # what is left is a combination of isotopes (including special
    # characters) and/or non-isotopes find things and remove them
    # order here matters because {'n', 'g', 'p', 't', 'd', 'a'} could all
    # be part of an isotope name or Webnucleo syntax name, in addition to
    # being valid themselves
    ret = []
    test_string = speciesString
    # start with proper nuclei first - [1:] because we don't want neutron
    for var in isotope_lut[1:]:
        num_var = test_string.count(var)
        for n in range(num_var):  # we might have more than 1 instance
            ret.append(var)
            test_string = test_string.replace(var, "")
    # check non_isotopes
    for var in rxn_to_WN_map:
        num_var = test_string.count(var)
        for n in range(num_var):  # we might have more than 1 instance
            ret.append(rxn_to_WN_map[var])
            test_string = test_string.replace(var, "")
    # this should be only special characters
    for var in specialCharZA:
        num_var = test_string.count(var)
        for n in range(num_var):  # we might have more than 1 instance
            ret.append(var)
            test_string = test_string.replace(var, "")
    # if there is anything left, then we didn't properly parse
    if test_string:
        errString = ("Didn't properly parse %s.  Ended up with %s"
                     % (speciesString, test_string))
        raise RuntimeError(errString)

    # # lone gammas
    # if speciesString.strip() == 'g':
    #     return 'gamma'
    # # what's left are either pure nuclei, pure special characters (e.g. 'a')
    # # or a mixture of nuclei and/or specials
    # # split them here - prune any empty strings that occur
    # split_species = [spec for spec in _specSplitter.split(speciesString)
    #                  if spec]
    # ret = []
    # for specString in split_species:
    #     # at this point, specString could still be a combination of
    #     # specials (e.g. 'an')
    #     # use the fact that get_Z_A will throw an KeyError if this is happens
    #     try:
    #         _ = get_Z_A(specString)
    #     except KeyError:
    #         # combination of specials; need to split them up, fix them,
    #         # and pass them back
    #         ret.extend([_fix_special_species(spec) for spec in specString])
    #     else:
    #         # this should be a valid entry in either Zdict or specialCharZA
    #         ret.append(_fix_special_species(specString))
    return PLUS.join(ret)


def form_rate_string(target, reactants, products, endState):
    """
    Takes the various parts of a reaction and forms the stochiometric version,
    i.e.
      target + reactants -> endState + products
    """
    lhs = PLUS.join(target + reactants)
    rhs = PLUS.join(endState + products)
    return lhs + TO + rhs
    # rstring = [target, ]
    # if reactants:
    #     rstring += [PLUS, reactants]
    # rstring.append(TO)
    # if endState:
    #     rstring.append(endState)
    # if products:
    #     rstring += [PLUS, products]

    # # this is of the form a + b -> c + d, with special character handling
    # rstring = ''.join(rstring)
    # return rstring


def _fix_special_species(species_string):
    """
    Fix the cases of special characters; i.e. t --> H3, a --> He4, etc.
    """
    if species_string in specialCharZA:
        Z, A = specialCharZA[species_string]
        # special handling of neutron
        if Z == 0:
            A = ''
        return "%s%s" % (element_lut[Z], A)

    return species_string
