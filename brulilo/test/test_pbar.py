from brulilo.util import progressbar
import time

nIsotopes = 50

mybar = progressbar.IntProgressBar("Getting Isotope data", nIsotopes)

for isotope in range(nIsotopes):
    time.sleep(0.1)
    mybar.update()
