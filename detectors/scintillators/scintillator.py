#This file is home to the parent Scintillator class
import numpy as np
import pandas as pd
from interfaces import IDetector, IParticle
from particles.photon import photon
import math
fano = 1 

class Scintillator(IDetector):
    def __init__(self):
        self.fano=fano
        
    def detects_batch(self, O, n, theta, phi, E, batch_size, gain):
        compton_bool = False
        compton = 0 #number of compton scattered photons. It will be used to recursively call detects_batch for multiple scatters.
       
        #Create a dataframe (2D array) of the values that will be obtained from the calculations for compton scatters
        df = pd.DataFrame(np.zeros((batch_size, 10)), columns=["compton_E", 
                                                               "photon_s_px", 
                                                               "photon_s_py", 
                                                               "photon_s_pz", 
                                                               "photon_s_E", 
                                                               "photon_s_theta", 
                                                               "photon_s_phi", 
                                                               "photon_s_x", 
                                                               "photon_s_y", 
                                                               "photon_s_z"])
        #Create a separate array for photoelectric electrons as it will have a different length
        photo_electrons = np.zeros(batch_size)
        photo_mask = np.zeros(batch_size).astype(bool)
        compton_mask = np.zeros(batch_size).astype(bool)
        positions = np.zeros(shape = (batch_size, 3))
        position_mask = np.zeros(batch_size).astype(bool)
        polar_positions = np.zeros(shape = (batch_size, 2)) #realistically this isn't needed since phi and theta are already inputs and this is just reassigning them 
        
        detected = 0
        #horrible array setups
        #print(n, self.a)
        a = np.tile(self.a, batch_size,)
        a = np.reshape(a, (batch_size, 3)) #reshaping the axis for vector maths (from shape (3,) to shape (batch_size,3))
        disc = np.sum((np.cross(n, a)* np.cross(n, a)), axis=1)*self.r**2 - (np.sum(((self.b - O)* np.cross(n, a)), axis=1)**2) #discriminant
        with np.errstate(divide='ignore', invalid='ignore'): #Distances towards cylinder surface intersection points
            d1 = (np.sum(np.cross(n, a)* np.cross((self.b-O), a), axis=1) - np.sqrt(disc))/(np.sum(np.cross(n, a)*np.cross(n, a), axis=1))
            d2 = (np.sum(np.cross(n, a)* np.cross((self.b-O), a), axis=1) + np.sqrt(disc))/(np.sum(np.cross(n, a)*np.cross(n, a), axis=1))
            d1 = np.reshape(d1, (batch_size,1))
            d1 = np.tile(d1, 3) #reshaping the distances ready for vector maths
            d2 = np.reshape(d2, (batch_size, 1))
            d2 = np.tile(d2, 3)
        
        dc1 = np.sum(a*(self.b-O), axis = 1)/np.sum(a*n, axis = 1) #distances towards cap intersections
        dc1 = np.reshape(dc1, (batch_size,1))
        dc1 = np.tile(dc1, 3)
        dc2 = np.sum(a*(self.t-O), axis = 1)/np.sum(a*n, axis = 1)
        dc2 = np.reshape(dc2, (batch_size,1))
        dc2 = np.tile(dc2, 3) 
        #batch_detected = np.sum(disc>=0) #flawed
        p1_cyl = O + d1*n #intersection with cylinder surface
        p2_cyl = O + d2*n
        angles = photon.gen_angles(E, batch_size)
        capcheck1 = (np.sum((O + n*dc1 - self.b)*(O + n*dc1 - self.b), axis=1)) <= self.r**2
        p1_cap = O + dc1*n
        p2_cap = O + dc2*n
        if E > 500:
            photo_csc = photon.photoelectric_csc_high(self.Z_n, E) #photoelectric crosssection for high energies
        else:
            photo_csc = photon.photoelectric_csc_mid(self.Z_n, E)
        compton_csc = photon.comptonscatter_csc(E, self.Z_n)
        total_csc = compton_csc + photo_csc
        threshold = photo_csc/total_csc
        #print(compton_csc, photo_csc)#FIXED THIS!! compton_csc was in mB so it was 1000x bigger than it's supposed to be!
        total = 0
        for i in range(batch_size):
            if (disc[i] >=0 and self.b[0] <= p1_cyl[i,0] <= self.t[0] and self.b[1] <= p1_cyl[i, 1] <= self.t[1] and self.b[2] <= p1_cyl[i,2] <= self.t[2]) or capcheck1[i] == True: #need to sort out  a check for the end caps
                detected += 1                #Going to compute interaction probability based off of photon cross sections
                if capcheck1[i] == False:
                    dist = abs(np.sum(d2[i,0] - d1[i,0])) #distance travelled inside of the detector (cm)
                    #These do not take into consideration whether or not the particle enters via the cap and leaves via the surface  & vise versa
                    p = p1_cyl[i]
                else:
                    dist = abs(np.sum(dc2[i, 0]-dc1[i, 0])) #Only want one component of them (theyre tiled for easier vector maths from earlier)
                    p = p1_cap[i]
                lifetime = 1/(self.n*total_csc*10**-14*3) #10^-14 since 1 barn  = 10^-24 cm^2, and c = 3*10^10 cm/s 
                #probability = 1 - np.e**(-dist/(lifetime*3*10**10))  #Need to check this since it seems very high?
                interaction_threshold = np.random.rand()
                interaction_type = np.random.rand() 
                #print(probability)
                #print(dist)
                if -np.log(interaction_threshold)*(lifetime*3*10**10)*10 <= dist: #this will find the distance the photon travels before interacting
                    #^needs fixing
                        #fixed?
                    #gonna have to redo this/check it? It's like a factor of 10 too small? Or at least 10 times smaller than previously? but dist is correct, so maybe it was just wrong before
                    total+=1
                    cos_phi = np.cos(phi[i])
                    dist_int = -np.log(interaction_threshold)*(lifetime*3*10**10)
                    #print(dist_int,"cm") This seems correct now
                    x_int = 10*dist_int * np.cos(theta[i]) * cos_phi + p[0]#The point of interaction 
                    y_int = 10*dist_int * np.sin(theta[i]) * cos_phi + p[1]
                    z_int = 10*dist_int * np.sin(phi[i]) + p[2] #This makes it kinda slow due to all of the trig
                    position = x_int, y_int, z_int
                    positions[i] = position
                    polar_position = theta[i], phi[i]
                    polar_positions[i] = polar_position
                    position_mask[i] = True
                    if interaction_type <= threshold: #my photopeak..... so so small
                        electron_E = (photon.photoelectric(theta[i], phi[i], E, x_int, y_int, z_int, dist_int, self.fano))
                        photo_electrons[i] = electron_E
                        photo_mask[i] = True
                        compton_bool = False
                    else:

                        #Get values from the comptonscatter function for the given photon
                        electron_energy, photon_px, photon_py, photon_pz, photon_energy, photon_theta, photon_phi = (photon.comptonscatter(theta[i], phi[i], E, x_int, y_int, z_int, dist_int, self.fano, angles))

                        #Insert them into the dataframe
                        df.loc[i] = electron_energy, photon_px, photon_py, photon_pz, photon_energy, photon_theta, photon_phi, x_int, y_int, z_int
                        compton_mask[i] = True

                        compton += 1
                        #the scattered photon from compton scattering
                        compton_bool = True

        #Apply masks to extract cases where interactions have occurred
        compton_indices = np.where(compton_mask)[0] # compton_df = df[compton_mask] 
        photo_electrons = photo_electrons[photo_mask]

        compton_electrons = np.zeros(compton)

        if compton > 0:
            #print(compton, "photons scattered in this batch") #Debug
            origins = df.loc[compton_indices, ["photon_s_x", "photon_s_y", "photon_s_z"]].to_numpy() # origins = compton_df[["photon_s_x", "photon_s_y", "photon_s_z"]].to_numpy()
            directions = df.loc[compton_indices, ["photon_s_px", "photon_s_py", "photon_s_pz"]].to_numpy() # directions = compton_df[["photon_s_px", "photon_s_py", "photon_s_pz"]].to_numpy()

            for i in range(compton):
                compton_bool = True
                electron_E = df.loc[compton_indices[i], "compton_E"] # electron_E = compton_df["compton_E"].iloc[i]
                j=0
                while compton_bool == True:
                    #print(electron_E, photon_s_E[i], j)
                    j +=1
                    electron = (self.detects_single(origins[i], directions[i], df.loc[compton_indices[i], "photon_s_theta"], df.loc[compton_indices[i], "photon_s_phi"], df.loc[compton_indices[i], "photon_s_E"], photo_electrons, d1[i], i, angles)) # electron = (self.detects_single(origins[i], directions[i], #recursively calling the method for each scattered photon. HORRIBLY inefficient but since the energies are different, each photon needs to be called inidividually
                    electron_E += electron[0]
                    compton_bool = electron[1]
                    if compton_bool == False:
                        break
                    directions[i] = (electron[2], electron[3], electron[4])/np.sqrt(electron[2]**2  + electron[3]**2 + electron[4]**2) #nomalise momence

                    df.loc[compton_indices[i], "photon_s_E"] = electron[5]      # compton_df["photon_s_E"].iloc[i] = electron[5]
                    df.loc[compton_indices[i], "photon_s_theta"] = electron[6]  # compton_df["photon_s_theta"].iloc[i] = electron[6]
                    df.loc[compton_indices[i], "photon_s_phi"] = electron[7]    # compton_df["photon_s_phi"].iloc[i] = electron[7]

                    #totals += electron[8]
                    #bad += electron [9]
                compton_electrons[i] = electron_E
                
        total_E_list = np.concatenate((photo_electrons, compton_electrons))
        scintillation_photons = total_E_list * self.ly
        positions = positions[position_mask]
        polar_positions = polar_positions[position_mask]
        return detected, scintillation_photons/100 * gain, total, total_E_list, positions, polar_positions

    #v This may end up being reduntant as I kind of want to use detects_batch again for all of the scatters within a batch, but I would need to change how getting angles works right now since it currently uses a single energy and not an array of energies. Would also need to make an array of crossections
    def detects_single(self, O, n, theta, phi, E, electron, t, i, angles): #this is just for compton photons   
        #Doesn't need to go through the detection algorithm, we already know it is in the detector (this will change when penetration into a detector is considered)
            #hi this is changing now
             #changing again for another detection algorithm
                 #the salami lid dont fit
                 #it works now but keep an eye on this.... I think too many internal scatters are happening
        disc = np.sum((np.cross(n, self.a)* np.cross(n, self.a)))*self.r**2 - (np.sum(((self.b - O)* np.cross(n, self.a)))**2) #discriminant
        
        with np.errstate(divide='ignore', invalid='ignore'): #Distance towards cylinder surface intersection points
            d1 = abs((np.sum(np.cross(n, self.a)* np.cross((self.b-O), self.a)) - np.sqrt(disc))/(np.sum(np.cross(n, self.a)*np.cross(n, self.a))))
            d2 = abs((np.sum(np.cross(n, self.a)* np.cross((self.b-O), self.a)) + np.sqrt(disc))/(np.sum(np.cross(n, self.a)*np.cross(n, self.a))))
        #should only really be d1 or d2 but they should be the same (disc ==0)
        #p1_cyl = O + n*d1
        #p2_cyl = O + n*d2
        #print(disc) #fmscl
        dc1 = abs(np.sum(self.a*(self.b-O))/np.sum(self.a*n)) #distance towards cap intersections
        dc2 = abs(np.sum(self.a*(self.t-O))/np.sum(self.a*n))
        #p1_cap = O + n*dc1
        #p2_cap = O + n*dc2
        if E > 500:
            photo_csc = photon.photoelectric_csc_high(self.Z_n, E) #photoelectric crosssection for high energies
        else:
            photo_csc = photon.photoelectric_csc_mid(self.Z_n, E)
        compton_csc = photon.comptonscatter_csc(E, self.Z_n)
        total_csc = compton_csc + photo_csc
        lifetime = 1/(self.n*total_csc*10**-14*3)
        threshold = photo_csc/total_csc #eventually this will approach 1 for low energy photons -> guaranteeing that comptonscatters will end with a photoelectri absorption (but also there's the probability of escaping which isn't considered here)
        random = np.random.rand()
        interaction_threshold = np.random.rand()
        if d1< d2:
            if d1 < dc1 or d1 <dc2:
                dist = d1
                #p = p1_cyl
            else:
                if dc1 < dc2:
                    dist = dc1
                    #p = p1_cap
                else:
                    dist = dc2
                    #p = p2_cap
        else:
            if d2 < dc1 or d2 <dc2:
                dist = d2
                #p = p2_cyl
            else:
                if dc1 < dc2:
                    dist = dc1
                    #p = p1_cap
                else:
                    dist = dc2
                    #p = p2_cap
        
        if -np.log(interaction_threshold)*(lifetime*3*10**10)*10 <= dist:#this is in mm (the distance is in cm*10)
            #print("dist", dist, "interaction distance", -np.log(1-interaction_threshold)*(lifetime*3*10**10)*10)
            #changed the interaction threshold calculation
            dist_int = -np.log(interaction_threshold)*(lifetime*3*10**10)
            cos_phi = np.cos(phi)
            x_int = 10*dist_int * np.cos(theta) * cos_phi + O[0]#The point of interaction in mm
            y_int = 10*dist_int * np.sin(theta) * cos_phi +O[1]
            z_int = 10*dist_int * np.sin(phi) + O[2]
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