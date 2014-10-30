from brulilo import Network
from brulilo import Reaction
from brulilo import Isotope

# CNO cycle
prot = Isotope("H", mass=1, ebin=0.0)
alfa = Isotope("He", mass=4, ebin=None)
C12 = Isotope("C", mass=12, ebin=None)
C13 = Isotope("C", mass=13, ebin=None)
N13 = Isotope("N", mass=13, ebin=None)
N14 = Isotope("N", mass=14, ebin=None)
N15 = Isotope("N", mass=15, ebin=None)
O15 = Isotope("O", mass=15, ebin=None)
O16 = Isotope("O", mass=16, ebin=None)
O17 = Isotope("O", mass=17, ebin=None)
F17 = Isotope("F", mass=17, ebin=None)
isotopes = [prot, alfa, C12, C13, N13, N14, N15,
            O15, O16, O17, F17]

c12pg = Reaction([(-1, C12), (1, prot), (1, N13)], None)
n13beta = Reaction([(-1, N13), (1, C13)], None)
c13pg = Reaction([(-1, C13), (1, prot), (1, N14)], None)
n14pg = Reaction([(-1, N14), (1, prot), (1, O15)], None)
o15beta = Reaction([(-1, O15), (1, N15)], None)
n15pg = Reaction([(-1, N15), (1, O16)], None)
n15pa = Reaction([(-1, N15), (-1, prot), (1, C12), (1, alfa)], None)
o16pg = Reaction([(-1, O16), (1, F17)], None)
f17beta = Reaction([(-1, F17), (1, O17)], None)
o17pa = Reaction([(-1, O17), (-1, prot), (1, N14), (1, alfa)], None)
rxns = [c12pg, n13beta, c13pg, n14pg, o15beta, n15pg,
        n15pa, o16pg, f17beta, o17pa]

net = Network(isotopes, rxns)
net.plot()
