import numpy as np

from particles.photon import photon
from detectors.naitl import NaITl

E = 662 #Energies are in keV, natural units are assumed, distances are in mm for the time being (for the sake of dimensions of detectors)

#Creating the detector once to avoid creating it 1000 times
detector = NaITl(100, 100, 0, 300, 300, 300)
detected = 0

n_photons = 10_000_000 # number of photons. also you can insert _ in a number without issue to break it up in python
batch_size = 100_000 # size per batch (i usually do ~10% of photons but its kind of up to you and your memory, at 10m i do 1% so 100k)

for i in range(0, n_photons, batch_size):
    theta = np.random.uniform(0, 2*np.pi, batch_size)
    phi = np.random.uniform(-np.pi/2, np.pi/2, batch_size)
    px = E * np.cos(theta) * np.cos(phi)
    py = E * np.sin(theta) * np.cos(phi)
    pz = E * np.sin(phi)
    directions = photon.batch_calculate_direction(pz, theta, phi)
    origins = np.zeros((batch_size, 3))  # assumes they all have the same origin of 0 0 0
    
    # DETECT THAT BATCH!!
    mask = detector.detects_batch(origins, directions)
    detected += np.sum(mask)

print(detected)
