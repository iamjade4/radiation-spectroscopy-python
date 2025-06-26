import numpy as np
from interfaces import IDetector, IParticle
from particles.photon import photon
import math
fano = 1 #This is intrinsic to scintillators in general, so will be defined here. Maybe if I make a parent class of Scintillator I could provide it there as the Fano factor tends to be 1 for all scintillators
Z_n = 53 #atomic number of Iodine. Can't remember if this is supposed to be for Iodine or for the activator. I assume Iodine but the result is so so smalllll
n = 1.474*10**22#number density
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
        self.Z_n = Z_n
        self.n = n

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
        tclose = np.reshape(tclose, (-1,1))
        tfar = np.reshape(tfar, (-1,1))
        positions = origins + np.multiply(tclose, directions)
        #print(np.shape(positions))
        #print(positions[0], tclose[0], directions[0])
        positions = positions.reshape(-1,3)
        tclose = tclose.flatten()
        tfar = tfar.flatten()
        x = np.where(tclose <= tfar, positions[:,0], None)
        y = np.where(tclose <= tfar, positions[:,1], None)
        z = np.where(tclose <= tfar, positions[:,2], None) #these are crazy big???
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
        total_csc = compton_csc + photo_csc
        threshold = photo_csc/total_csc
        #print(compton_csc, photo_csc) the photoelectric crosssection seems sooooo small for high energies but i guess that makes sense
        total = 0
        for i in range(batch_size):
            if tclose[i] <= tfar[i]:
                #Going to compute interaction probability based off of photon cross sections
                dist = (tfar[i]-tclose[i])/10 #distance travelled inside of the detector (cm)
                lifetime = 1/(n*total_csc*10**-14*3) #10^-14 since 1 barn  = 10^-24 cm^2, and c = 3*10^10 cm/s 
                #probability = 1 - np.e**(-dist/(lifetime*3*10**10))  #Need to check this since it seems very high?
                interaction_threshold = np.random.rand()
                interaction_type = np.random.rand() 
                #print(probability)
                if -np.log(1-interaction_threshold)*(lifetime*3*10**10) <= dist: #this will find the distance the photon travels before interacting
                    total+=1
                    cos_phi = np.cos(phi[i])
                    dist_int = -np.log(1-interaction_threshold)*(lifetime*3*10**10)
                    x_int = dist_int * np.cos(theta[i]) * cos_phi + x[i]#The point of interaction 
                    y_int = dist_int * np.sin(theta[i]) * cos_phi + y[i]
                    z_int = dist_int * np.sin(phi[i]) + z[i] #This makes it kinda slow due to all of the trig
                    if interaction_type <= threshold: #my photopeak..... so so small
                        electron_E = (photon.photoelectric(thetas[i], phis[i], E, x_int, y_int, z_int, tclose_det[i], self.fano))  
                        electrons.append(electron_E)
                        compton_bool = False
                    else:
                        electron = (photon.comptonscatter(thetas[i], phis[i], E, x_int, y_int, z_int, tclose_det[i], self.fano, angles))
                        compton_E.append(electron[0])
                        compton += 1
                        photon_s_px.append(electron[1])
                        photon_s_py.append(electron[2])
                        photon_s_pz.append(electron[3])
                        photon_s_E.append(electron[4])
                        photon_s_theta.append(electron[5])
                        photon_s_phi.append(electron[6])
                        #totals += electron[7]
                        #bad += electron[8] #debugging
                        photon_s_x.append(x_int)
                        photon_s_y.append(y_int)
                        photon_s_z.append(z_int) #this sucks
                        #the scattered photon from compton scattering
                        compton_bool = True
        if compton > 0:
            #print(compton, "photons scattered in this batch") #Debug
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
                    #totals += electron[8]
                    #bad += electron [9]
                electrons.append(electron_E)
        return tclose <= tfar, electrons, total#The electron is returning a 2d array of ALL of the different electron energies. For now it will all be 662        
            
    def detects_single(self, origin, direction, theta, phi, E, electron, t, i, angles): #this is just for compton photons        
        #Doesn't need to go through the detection algorithm, we already know it is in the detector (this will change when penetration into a detector is considered)
            #hi this is changing now
        xc, yc, zc, xf, yf, zf = self._bounds
        bound_min = np.array([xc, yc, zc])
        bound_max = np.array([xf, yf, zf])

        with np.errstate(divide='ignore', invalid='ignore'):
            t_min = (bound_min - origin) / direction
            t_max = (bound_max - origin) / direction
        t_min = np.where(direction == 0, -np.inf, t_min)
        t_max = np.where(direction == 0, np.inf, t_max)
        tclose = np.max(t_min)
        tfar = np.min(t_max)

        position = origin + tclose*direction #This needs redoing
        #print(positions, tclose,tfar, directions)
        x = position[0]
        y = position[1]
        z = position[2]
        #x = np.where(tclose <= tfar, positions[:,0], None)
        #y = np.where(tclose <= tfar, positions[:,1], None)
        #z = np.where(tclose <= tfar, positions[:,2], None)
        dist = math.sqrt((x-origin[0])**2 + (y-origin[1])**2 + (z-origin[2])**2)/10  #YEAHHHH
        
        #print(dist)
        if E > 500:
            photo_csc = photon.photoelectric_csc_high(Z_n, E) #photoelectric crosssection for high energies
        else:
            photo_csc = photon.photoelectric_csc_mid(Z_n, E)
        compton_csc = photon.comptonscatter_csc(E)
        total_csc = compton_csc + photo_csc
        lifetime = 1/(n*total_csc*10**-14*3)
        threshold = photo_csc/total_csc #eventually this will approach 1 for low energy photons -> guaranteeing that comptonscatters will end with a photoelectri absorption (but also there's the probability of escaping which isn't considered here)
        random = np.random.rand()
        interaction_threshold = np.random.rand()
        if -np.log(1-interaction_threshold)*(lifetime*3*10**10) <= dist:
            dist_int = -np.log(1-interaction_threshold)*(lifetime*3*10**10)
            cos_phi = np.cos(phi)
            x_int = dist_int * np.cos(theta) * cos_phi + x#The point of interaction 
            y_int = dist_int * np.sin(theta) * cos_phi + y
            z_int = dist_int * np.sin(phi) + z
            if random <= threshold: #my photopeak..... so so small
                electron_E = (photon.photoelectric(theta, phi, E, x_int, y_int, z_int, t, self.fano))  
                compton_bool = False
                return electron_E, compton_bool
            else:
                electron = (photon.comptonscatter(theta, phi, E, x_int, y_int, z_int, t, self.fano, angles))
                compton_bool = True
                electron_E = electron[0]    
                return electron_E, compton_bool, electron[1], electron[2], electron[3], electron[4], electron[5], electron[6]
        else:
            return(0, 0)