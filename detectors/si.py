import numpy as np
from interfaces import IDetector, IParticle
from particles.photon import photon
electrons = []
fano = 0.129 #General value for silico detectors. It does vary with whether or not it is cryogenically cooled though so this is just a general estimate

class Si(IDetector):
    def __init__(self, height, width, thickness, n, position):
        #Lowercase are the minumum bound position and uppercase are the dimensions
        #So far this only allows for cuboidal detectors, I will have to redefine it for more complex geometry
        #This is a silicon slab detector so a cuboidal geometry actually suits it here
        self.h = height
        self.w = width
        self.d = thickness/1000
        self.n = n
        self.o = position
        self._bounds = self._calculate_bounds()
        self.fano = fano #This allows multiple detector types to be called at once by removing the need for fano to remain a single value, it is now inherent to the detector's instance

    def _calculate_bounds(self):
        #thickness is in x axis by default (theta = 0, phi = 0)
        self.x = self.o[0]
        self.y = self.o[1] - self.w/2
        self.z = self.o[2] - self.h/2
        self.X = self.o[0] + self.d
        self.Y = self.o[1] + self.w/2
        self.Z = self.o[2] + self.h/2
        if self.n[1] != 0:
            self.theta  = np.arccos(self.n[0]/self.n[1])
        else:
            self.theta = 0
        self.phi = np.arcsin(self.n[2])
        sin_theta = np.sin(self.theta)
        cos_theta = np.cos(self.theta)
        sin_phi = np.sin(self.phi)
        cos_phi = np.cos(self.phi)
        Rz = np.array([cos_theta, -sin_theta, 0,
                      sin_theta, cos_theta, 0,
                      0, 0, 1])
        Rz = np.reshape(Rz, (3,3))
        Ry = np.array([cos_phi, 0, -sin_phi,
                       0, 1, 0,
                       sin_phi, 0, cos_phi])
        Ry = np.reshape(Ry, (3,3))
        self.x, self.y, self.z = np.dot(Rz,np.dot(Ry,[self.x, self.y, self.z]))
        self.X, self.Y, self.Z = np.dot(Rz,np.dot(Ry,[self.X, self.Y, self.Z]))
        
        
        return (
            self.x, self.y, self.z,
            self.X,
            self.Y,
            self.Z
        )

        
    def detects_batch(self, origins, directions, theta, phi, E, batch_size):
        electrons = []
        xc, yc, zc, xf, yf, zf = self._bounds
        bound_min = np.array([xc, yc, zc])
        bound_max = np.array([xf, yf, zf])
        detected = 0
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
        for i in range(batch_size):
            if tclose[i] <= tfar[i]:
                detected +=1
                electron_E = (photon.photoelectric(thetas[i], phis[i], E, x[i], y[i], z[i], tclose_det[i], self.fano))  
                electrons.append(round(electron_E)) #rounding to an int here. Realistically, it will just be put into a channel number for spectroscopy but rounding is easiest for now
                #Silicon detectors are generally more suited to direct detection of charged particles, and not necessarily of photons due to their much smaller detector volume compared to that of scintillators. For now however, for the sake of comparison of energy resolutions, I will keep the photoelectric method here
        
        #print(positions)
        return detected, electrons, detected #The electron is returning a 2d array of ALL of the different electron energies. For now it will all be 662        
            
        