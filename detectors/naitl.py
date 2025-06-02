import numpy as np
from interfaces import IDetector, IParticle
from particles.photon import photon
fano = 1 #This is intrinsic to scintillators in general, so will be defined here. Maybe if I make a parent class of Scintillator I could provide it there as the Fano factor tends to be 1 for all scintillators
Z_n = 53 #atomic number of Iodine. Can't remember if this is supposed to be for Iodine or for the activator. I assume Iodine but the result is so so smalllll
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

        
    def detects_batch(self, origins, directions, theta, phi, E, batch_size):
        compton_bool = False
        compton = 0 #number of compton scattered photons. It will be used to recursively call detects_batch for multiple scatters.
        photon_s_px = []
        photon_s_py = []
        photon_s_pz = []
        photon_s_E = []
        photon_s_theta = []
        photon_s_phi = []
        photon_s_x = []
        photon_s_y = []
        photon_s_z = []
        electrons = []
        compton_E = []
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
        positions = origins + np.dot(tclose, directions)
        positions = positions.reshape(-1,3)
        
        x = np.where(tclose <= tfar, positions[:,0], None)
        y = np.where(tclose <= tfar, positions[:,1], None)
        z = np.where(tclose <= tfar, positions[:,2], None)
        thetas = np.where(tclose <= tfar, theta, None)
        phis = np.where(tclose <= tfar, phi, None)
        batch_detected = np.sum(tclose <= tfar)
        tclose_det = np.where(tclose <= tfar, tclose, None)
        angles = photon.gen_angles(E, batch_size)
        if E > 500:
            photo_csc = photon.photoelectric_csc_high(Z_n, E) #photoelectric crosssection for high energies
        else:
            photo_csc = photon.photoelectric_csc_mid(Z_n, E)
        compton_csc = photon.comptonscatter_csc(E)
        total = compton_csc + photo_csc
        threshold = photo_csc/total
        #print(compton_csc, photo_csc) the photoelectric crosssection seems sooooo small for high energies but i guess that makes sense
        for i in range(batch_size):
            if tclose[i] <= tfar[i]:
                #Going to compute interaction probability based off of photon cross sections
                random = np.random.rand()
                if random <= threshold: #my photopeak..... so so small
                    electron_E = (photon.photoelectric(thetas[i], phis[i], E, x[i], y[i], z[i], tclose_det[i], self.fano))  
                    electrons.append(electron_E)
                    compton_bool = False
                else:
                    electron = (photon.comptonscatter(thetas[i], phis[i], E, x[i], y[i], z[i], tclose_det[i], self.fano, angles))
                    compton_E.append(electron[0])
                    compton += 1
                    photon_s_px.append(electron[1])
                    photon_s_py.append(electron[2])
                    photon_s_pz.append(electron[3])
                    photon_s_E.append(electron[4])
                    photon_s_theta.append(electron[5])
                    photon_s_phi.append(electron[6])
                    photon_s_x.append(x[i])
                    photon_s_y.append(y[i])
                    photon_s_z.append(z[i]) #this sucks
                    #the scattered photon from compton scattering
                    compton_bool = True
        if compton > 0:
            print(compton)
            origins = photon_s_x, photon_s_y, photon_s_z
            origins = np.transpose(origins)
            directions = photon_s_px, photon_s_py, photon_s_pz
            directions = np.transpose(directions)
            for i in range(compton):
                compton_bool = True
                electron_E = compton_E[i]
                j=0
                while compton_bool == True:
                    #print(electron_E, photon_s_E[i], j)
                    j +=1
                    electron = (self.detects_single(origins[i], directions[i], photon_s_theta[i], photon_s_phi[i], photon_s_E[i], electrons, tclose_det[i], i, angles)) #recursively calling the method for each scattered photon. HORRIBLY inefficient but since the energies are different, each photon needs to be called inidividually
                    electron_E += electron[0]
                    compton_bool = electron[1]
                    if compton_bool == False:
                        break
                    directions[i] = electron[2], electron[3], electron[4]
                    photon_s_E[i]= electron[5]
                    photon_s_theta[i] = electron[6]
                    photon_s_phi[i] = electron[7]
                    #print(electron_E + photon_s_E[i]) #this should always be 662. But it isnt!
                electrons.append(electron_E)
        return tclose <= tfar, electrons #The electron is returning a 2d array of ALL of the different electron energies. For now it will all be 662        
            
    def detects_single(self, origins, directions, thetas, phis, E, electron, t, i, angles): #this is just for compton photons        
        #Doesn't need to go through the detection algorithm, we already know it is in the detector (this will change when penetration into a detector is considered)
        if E > 500:
            photo_csc = photon.photoelectric_csc_high(Z_n, E) #photoelectric crosssection for high energies
        else:
            photo_csc = photon.photoelectric_csc_mid(Z_n, E)
        compton_csc = photon.comptonscatter_csc(E)
        total = compton_csc + photo_csc
        threshold = photo_csc/total #eventually this will approach 1 for low energy photons -> guaranteeing that comptonscatters will end with a photoelectri absorption (but also there's the probability of escaping which isn't considered here)
        random = np.random.rand()
        if random <= threshold: #my photopeak..... so so small
            electron_E = (photon.photoelectric(thetas, phis, E, origins[0], origins[1], origins[2], t, self.fano))  
            compton_bool = False
            return electron_E, compton_bool
        else:
            electron = (photon.comptonscatter(thetas, phis, E, origins[0], origins[1], origins[2], t, self.fano, angles))
            compton_bool = True
            electron_E = electron[0]    
            return electron_E, compton_bool, electron[1], electron[2], electron[3], electron[4], electron[5], electron[6]