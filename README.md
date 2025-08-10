# Radiation Spectroscopy (Python)
WIP project that simulates particle detection and radiation spectroscopy. It will, at some point:
- Accurately simulate real world effects that decrease detector resolution (i.e. bremsstrahlung loss, quantum efficiency statistics etc.)
- Simulate various gamma-ray interactions (compton scattering, pair production, photoelectric effect)

Once this proof of concept is sufficient, the plan is to port it to a 3D interactive realtime environment in Godot.

It uses an isotropic source with particle tracking. Gamma ray sources will produce particles with the class "Photon" and will carry px py pz and x y z information. They may also carry energy information as an explicit variable to save needing to convert momentum to energy, as well as a time signature. The source could create an array element (a particle) every 1/A seconds, where A is the source activity in bq. Hopefully it should be able to handle activities up to 1kbq, and ideally should handle higher activities but I am quite skeptical about this.

I will start by using an open source that is immune to self-absorption, and acts as a pure gamma source. Along the line I would like to take into consideration the effects of self-absorption, especially for charged particle sources where this may be significant. 

## Installation
This program is primarily developed and tested on Linux. The installation process will look a little different depending on your OS. If you run on a different OS you may experience previously unseen bugs. By all means, report these as Issues if you come across them. 

Download python **3.13.5** or later (earlier might work, but this is the version currently developed on). If you are on windows make sure to add python to your PATH.

Download the following python dependencies (either through pip or a package manager):
- Matplotlib
- Numpy
- Scipy
- PyQT5 (PyQT5 is used over 6 due to matplotlib compatibility)
- pandas

Currently we do not do releases, either download this repository as a .zip or `git clone https://github.com/iamjade4/radiation-spectroscopy-python.git` into your folder of choice. 

Run main.py, you may need to do so in the console by typing `python main.py` in the project's folder. 

## Contribution
Follow the Installation Instructions above. Fork the repository, make your changes, pull request.

Submit Issues if you come across them, even if they are OS specific. If you know that it is OS specific please leave your OS in the issue. Do not submit Issues for being unable to launch the program due to you lacking a certain dependency.