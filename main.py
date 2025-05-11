import numpy as np

from particles.photon import photon
from detectors.naitl import NaITl
from detectors.si import Si
from particles.electron import electron
import matplotlib.pyplot as plt

E = 662 #Energies are in keV, natural units are assumed, distances are in mm for the time being (for the sake of dimensions of detectors)
#1eV = 1.66e10-19J

#Creating the detector once to avoid creating it 1000 times
energies = []

detector1 = NaITl(100, 100, 0, 60, 60, 60) #Recently changed this to be 6x6x6 cm to be more realistic. Expect much smaller detection counts
detector2 = Si(100, 0, 100, 60, 60, 2) #Silicon detectors tend to be THIN (hence the 2mm depth) 
detected1 = 0
detected2 = 0

n_photons = 10_000_000 # number of photons. also you can insert _ in a number without issue to break it up in python
batch_size = 100_000 # size per batch (i usually do ~10% of photons but its kind of up to you and your memory, at 10m i do 1% so 100k)

for i in range(0, n_photons, batch_size):
    batch_detected1 = 0
    batch_detected2 = 0
    theta = np.random.uniform(0, 2*np.pi, batch_size)
    phi = np.random.uniform(-np.pi/2, np.pi/2, batch_size)
    px = E * np.cos(theta) * np.cos(phi)
    py = E * np.sin(theta) * np.cos(phi)
    pz = E * np.sin(phi)
    directions = photon.batch_calculate_direction(theta, phi)
    origins = np.zeros((batch_size, 3))  # assumes they all have the same origin of 0 0 0
    
    # DETECT THAT BATCH!!
    returns = [detector1.detects_batch(origins, directions, theta, phi, E), detector2.detects_batch(origins, directions, theta, phi, E)]
    
    mask1 = returns[0][0] #This is the bulk of bool (BULK OF BOOL!!!) ie number of detections in this batch
    detected1 += np.sum(mask1)
    batch_detected1 = np.sum(mask1)
    
    mask2 = returns[1][0]
    detected2 += np.sum(mask2)
    batch_detected2 = np.sum(mask2)
energies = np.concatenate((returns[0][1][:], returns[1][1][:]))#This is giving energy two different arrays, one for each detector. I will combine them just to superimpose the spectra into a single output, but realistically they should remain separate and have separate outputs
           #Concatenate combines the two arrays                         
print(detected1, detected2)

fig, ax = plt.subplots()
ax.hist(energies, bins=1024, range=[0,1000], histtype='step')
plt.show()