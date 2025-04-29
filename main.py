import numpy as np
import random as rd

from particles.photon import photon
from detectors.naitl import NaITl


photons = []
detected = 0
E = 662 #Energies are in keV, natural units are assumed, distances are in mm for the time being (for the sake of dimensions of detectors)

#Creating the detector once to avoid creating it 1000 times
detector = NaITl(100, 100, 0, 300, 300, 300)

for i in range(100000):
    #Randomly generated angles however, it will be generating a duplicate for 0 & 2pi, I could find a solution for this (go up to 2pi-1x10^-16 for flt16)
    theta = rd.uniform(0, 2*np.pi)
    phi = rd.uniform(-np.pi/2, np.pi/2)
    px = E*np.cos(theta)*np.cos(phi)
    py = E*np.sin(theta)*np.cos(phi)
    pz = E*np.sin(phi)
    #General momentum calculation, energy is based off of Cs 137 gammas. An actual source will not be mono-energetic (other radiation branches)
    Photon = photon(px, py, pz, theta, phi, 0, 0, 0, 0)
    #position at 0 0 0 since for now all photons are originating at the origin of the source
    photons.append(Photon)
    #Need to figure out sleep command so that it actually generates particles in a set timeframe. Also, the particles shouldn't be equally spaced out time-wise, they're supposed to represent a randomly disintegrating source
    if detector.detects(Photon):
        detected += 1

print(detected)
