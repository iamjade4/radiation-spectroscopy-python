import numpy as np
from interfaces import IDetector, IParticle
from particles.photon import photon
electrons = []
fano = 1 #This is intrinsic to scintillators in general, so will be defined here. Maybe if I make a parent class of Scintillator I could provide it there as the Fano factor tends to be 1 for all scintillators

class NaITl(IDetector):
    def __init__(self, x, y, z, X, Y, Z):
        #Lowercase are the minumum bound position and uppercase are the dimensions
        #So far this only allows for cuboidal detectors, I will have to redefine it for more complex geometry
        self.x = x
        self.y = y
        self.z = z
        self.X = X
        self.Y = Y
        self.Z = Z
        self._bounds = self._calculate_bounds()
        self.fano = fano #This allows multiple detector types to be called at once by removing the need for fano to remain a single value, it is now inherent to the detector's instance

    def _calculate_bounds(self):
        return (
            self.x, self.y, self.z,
            self.x + self.X,
            self.y + self.Y,
            self.z + self.Z
        )

        
    def detects_batch(self, origins, directions, theta, phi, E):
        xc, yc, zc, xf, yf, zf = self._bounds
        bound_min = np.array([xc, yc, zc])
        bound_max = np.array([xf, yf, zf])

        with np.errstate(divide='ignore', invalid='ignore'):
            t_min = (bound_min - origins) / directions
            t_max = (bound_max - origins) / directions

        t_min = np.where(directions == 0, -np.inf, t_min)
        t_max = np.where(directions == 0, np.inf, t_max)
        tclose = np.max(t_min, axis=1)
        tfar = np.min(t_max, axis=1)
        #Origins is (100000, 3) and I want to pull (100000,x) out of it
        positions = origins + np.dot(tclose, directions)
        positions = positions.reshape(-1,3)
        
        x = np.where(tclose <= tfar, positions[:,0], None)
        y = np.where(tclose <= tfar, positions[:,1], None)
        z = np.where(tclose <= tfar, positions[:,2], None)
        thetas = np.where(tclose <= tfar, theta, None)
        phis = np.where(tclose <= tfar, phi, None)
        batch_detected = np.sum(tclose <= tfar)
        tclose_det = np.where(tclose <= tfar, tclose, None)
        for i in range(100000):
            if tclose[i] <= tfar[i]:
                random = np.random.rand()
                if random <= 0.5:
                    electron_E = (photon.photoelectric(thetas[i], phis[i], E, x[i], y[i], z[i], tclose_det[i], self.fano))  
                else:
                    electron_E = (photon.comptonscatter(thetas[i], phis[i], E, x[i], y[i], z[i], tclose_det[i], self.fano)) #This makes it a 50/50 between making a photoelectron and a compton electron. This is not physical because I Am Tired
                electrons.append(round(electron_E)) #rounding to an int here. Realistically, it will just be put into a channel number for spectroscopy but rounding is easiest for now

        #print(positions)
        
        return tclose <= tfar, electrons #The electron is returning a 2d array of ALL of the different electron energies. For now it will all be 662        
            
        