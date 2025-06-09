from particles.photon import photon
from detectors.naitl import NaITl
from detectors.si import Si
from qtwidget import MyWidget
#from particles.electron import electron
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from PyQt5 import QtCore as Qtcre, QtWidgets as Qtwdgt, QtGui as Qtgui
import sys

SIMULATION_CONFIG = {
    "E": 662,  # Energies are in keV, natural units are assumed, distances are in mm for now (for the sake of dimensions of detectors), 1eV = 1.66e10-19J
    "n_photons": 10_000_000,
    "batch_size": 100_000,
    "bins": 1024,
    "hist_range": (0, 1024)
}
plt.rcParams['axes.xmargin'] = 0
totals = 0
bad = 0 #Debug variables
def simulate(n_photons: int, batch_size: int, E: float, detectors: list, totals, bad):
    total = [0] * len(detectors)
    detected_counts = [0] * len(detectors)
    energies = [[] for i in detectors]

    for i in range(0, n_photons, batch_size):
        theta = np.random.uniform(0, 2 * np.pi, batch_size)
        phi = np.random.uniform(-np.pi/2, np.pi/2, batch_size)
        directions = photon.batch_calculate_direction(theta, phi)
        origins = np.zeros((batch_size, 3))

        for idx, det in enumerate(detectors):
            mask, batch_energies, totals, bad, total_dtd = det.detects_batch(origins, directions, theta, phi, E, batch_size, totals, bad)
            total[idx] +=np.sum(total_dtd)
            detected_counts[idx] += np.sum(mask)
            energies[idx].extend(batch_energies) # this adds the cumulative data, previously you were only plotting the last batch (unless thats what you wanted to do in which case, sorry)
        print(i/batch_size+1,"/",(n_photons/batch_size), "batches simulated", end='\r')
    print(bad/totals *100, "Percentage of 'bad' scatters")
    return detected_counts, energies, total

def plot_spectra(energies: list, bins: int = 1024, energy_range=(0, 1000)):
    num_det = len(energies)
    fig, axs = plt.subplots(num_det, 1, figsize=(10, 4 * num_det), squeeze=False)
    for i, (ax, detector_energies) in enumerate(zip(axs.flat, energies)):
        ax.hist(detector_energies, bins=bins, range=energy_range, histtype='step')
    plt.tight_layout()
    plt.show() 
    return fig


def main():
    cfg = SIMULATION_CONFIG
    detectors = [
        NaITl(100, 100, 0, 60, 60, 60),  # 6x6x6cm NaITl
        Si(100, 0, 100, 60, 60, 2)       # 6x6x0.2cm Si (silicon detector is thin)
    ]
    detected_counts, energies, total = simulate(cfg["n_photons"], cfg["batch_size"], cfg["E"], detectors, totals, bad)
    for idx, count in enumerate(detected_counts, start=1):
        print(total[idx-1]/count *100, "% Detector efficiency")
        print(f"detected {count} photons")
    fig = plot_spectra(energies, bins=cfg["bins"], energy_range=cfg["hist_range"])
    if not Qtwdgt.QApplication.instance():
        app = Qtwdgt.QApplication([])
    else:
        app = Qtwdgt.QApplication.instance()
    widget = MyWidget(fig)
    widget.resize(1400, 800)
    widget.show()
    
    exit_code = app.exec()
    sys.exit(exit_code)
        

if __name__ == '__main__':
    main()
