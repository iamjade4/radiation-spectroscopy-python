This is a WIP project to simulate particle detection and radiation spectroscopy. I would like to eventually port this into a 3D interactive realtime environment using the Godot engine, but for now I am making crude mockups in python as a proof of concept. I would like to simulate the effects that decrease the resolution in detectors (i.e bremsstrahlung loss, quantum efficiency statistics etc.) as well as simulate various gamma-ray interactions (compton scattering, pair production, photoelectric effect).

Isotropic source with particle tracking. Gamma ray sources will produce particles with the class "Photon" and will carry px py pz and x y z information. They may also carry energy information as an explicit variable to save needing to convert momentum to energy, as well as a time signature. The source could create an array element (a particle) every 1/A seconds, where A is the source activity in bq. Hopefully it should be able to handle activities up to 1kbq, and ideally should handle higher activities but I am quite skeptical about this.

I will start by using an open source that is immune to self-absorption, and acts as a pure gamma source. Along the line I would like to take into consideration the effects of self-absorption, especially for charged particle sources where this may be significant. 

Required python dependencies:
- Matplotlib
- Numpy
- Scipy
- PyQT5