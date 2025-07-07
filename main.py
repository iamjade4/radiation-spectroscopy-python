#GUI
import sys
from PyQt5.QtCore import QSize, Qt
from PyQt5.QtWidgets import QApplication, QWidget, QMainWindow, QPushButton, QComboBox, QGridLayout, QSlider, QLabel, QProgressBar
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as Figcan
from particles.photon import photon
from detectors.naitl import NaITl
from detectors.si import Si
import matplotlib
import matplotlib.pyplot as plt
import numpy as np


plt.rcParams['axes.xmargin'] = 0

class MainWindow(QMainWindow):
    def __init__(self):#
        #Defaults
        self.E = 662 #keV
        self.batch_size =100_000
        self.n_photons = 10_000_000
        super().__init__()
        self.setWindowTitle("Radiation Spectroscopy Sim")
        layout = QGridLayout()
        self.layout= layout
        self.running = False

        #First box and text
        photon_drop = QComboBox(self)
        photon_drop.addItem("1_000_000")
        photon_drop.addItem("10_000_000")
        photon_drop.addItem("100_000_000")
        photon_drop.currentTextChanged.connect(self.n_photons_drop)
        photon_drop.setFixedHeight(30)
        photon_drop.setCurrentText("10_000_000")
        photons_label = QLabel("Number of photons (Default = 10,000,000):")
        photons_label.setFixedHeight(30)
        layout.addWidget(photons_label,0,0)
        layout.addWidget(photon_drop,1,0)
        
        #Second box and text
        batch_select = QComboBox(self)
        batch_select.addItem("1_000")
        batch_select.addItem("10_000")
        batch_select.addItem("100_000")
        batch_select.currentTextChanged.connect(self.batch)
        batch_select.setFixedHeight(30)
        batch_select.setCurrentText("100_000")
        batches_label = QLabel("Number of photons per batch (Default = 100,000):")
        batches_label.setFixedHeight(30)
        layout.addWidget(batches_label,2,0)
        layout.addWidget(batch_select,3,0)
        
        #Energy slider
        self.energy = QSlider(Qt.Orientation.Horizontal, self)
        self.energy.setRange(1,1000) # changed the lowest value to 1 because 0 crashes the program -K
        self.energy.setSingleStep(1)
        self.energy.setFixedHeight(30)
        self.energy.setValue(662)
        self.E = self.energy.value()
        self.energy.valueChanged.connect(self.energy_update)
        layout.addWidget(self.energy,5,0)
        
        #Energy label
        self.E_box = QLabel("Energy: "+str(self.E)+"keV (Default = 662keV)")
        self.E_box.setFixedHeight(30)
        layout.addWidget(self.E_box,4,0)
        
        #Start button
        start_button = QPushButton("Start")
        start_button.clicked.connect(self.main) #This runs the main that was from the backend
        start_button.setFixedHeight(30)
        layout.addWidget(start_button,6,0)
        #layout.addStretch(1)
        layout.setRowStretch(layout.rowCount(), 1)
        layout.setColumnStretch(layout.columnCount(), 1)
        ### creating widget with layout ###
        widget = QWidget()
        widget.setLayout(layout)
        self.setCentralWidget(widget)

    def n_photons_drop(self, n):
        self.n_photons = int(n)
        print("photons:", self.n_photons)

    
    def display_fig(self, fig):
        figure = Figcan(fig)
        height = self.layout.rowCount()
        self.layout.addWidget(figure,0,1,height,3)
        
    def batch(self, b):
        self.batch_size = int(b)
        print("batch size:", self.batch_size)
    
    def energy_update(self, E):
        self.E_box.setText("Energy: "+str(E)+"keV (Default = 662keV)")
        self.E = int(E)        
        
    def simulate(self, n_photons: int, batch_size: int, E: float, detectors: list):
        if E==0:
            print("WARNING. 0 PASSED INTO E.")

        total = [0] * len(detectors)
        detected_counts = [0] * len(detectors)
        energies = [[] for i in detectors]

        for i in range(0, n_photons, batch_size):
            theta = np.random.uniform(0, 2 * np.pi, batch_size)
            phi = np.random.uniform(-np.pi/2, np.pi/2, batch_size)
            directions = photon.batch_calculate_direction(theta, phi)
            origins = np.zeros((batch_size, 3))

            for idx, det in enumerate(detectors):
                mask, batch_energies, total_dtd = det.detects_batch(origins, directions, theta, phi, E, batch_size)
                total[idx] +=np.sum(total_dtd)
                detected_counts[idx] += np.sum(mask)
                energies[idx].extend(batch_energies) # this adds the cumulative data, previously you were only plotting the last batch (unless thats what you wanted to do in which case, sorry)
            print(i/batch_size+1,"/",(n_photons/batch_size), "batches simulated", end='\r')
            progress=  (i/batch_size+1)/(n_photons/batch_size)
            app.processEvents()
            self.batch_prog.setValue(int(progress*100))
        return detected_counts, energies, total

    def plot_spectra(self, energies: list, bins: int = 1024, energy_range=(0, 1000)):
        num_det = len(energies)
        fig, axs = plt.subplots(num_det, 1, figsize=(10, 4 * num_det), squeeze=False)
        for i, (ax, detector_energies) in enumerate(zip(axs.flat, energies)):
            ax.hist(detector_energies, bins=bins, range=energy_range, histtype='step')
        plt.tight_layout()
        #plt.show() 
        return fig


    def main(self):
        if self.running == False:
            self.patience = 1
            self.running = True
            E = self.E
            batch_size = self.batch_size
            n_photons = self.n_photons
            self.batch_prog = QProgressBar()
            self.layout.addWidget(self.batch_prog,7,0)
            self.layout.setRowStretch(self.layout.rowCount(), 1) #This is what resizes the graph when ran again but it stops the left column from restretching
            detectors = [
                NaITl(100, 100, 0, 60, 60, 60),  # 6x6x6cm NaITl
                Si(100, 0, 100, 60, 60, 2)       # 6x6x0.2cm Si (silicon detector is thin)
            ]
            detected_counts, energies, total = self.simulate(n_photons, batch_size, E, detectors)
            for idx, count in enumerate(detected_counts, start=1):
                print(total[idx-1]/count *100, "% Detector efficiency")
                print(f"detected {count} photons")
            fig = self.plot_spectra(energies, bins=1024, energy_range=(0, 1024))#removed calling config so that it can be easier edited at this point by the user
            self.display_fig(fig)
            self.running = False
            if self.patience == 0: #Remove the label if it has been made
                self.layout.removeWidget(self.patient)
                self.patient.deleteLater()
                self.patient = None
        else:
            if self.patience != 0:
                self.patience = 0
                self.patient = QLabel("<a href=\"https://jokerquotes.xyz\">'Be Patient'</a>" ) #Happens if a user presses "Start" while the program is currently running
                self.patient.setOpenExternalLinks(True)
                self.layout.addWidget(self.patient, 8, 0)
                
    def closeEvent(self,event):
        raise SystemExit(0)

app = QApplication(sys.argv)
window = MainWindow()
window.resize(1000,800)
window.show()
app.exec()  