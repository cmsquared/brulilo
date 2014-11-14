import brulilo
from brulilo import Network
import os.path

# this is a weak connection now, should probably generalize
rxn_file = os.path.join(os.path.dirname(brulilo.__file__),
                        'test/testRxns.txt')

net = Network.from_rxn_file(rxn_file)
#net.plot()

net.pprint()
