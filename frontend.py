import sys
from multiprocessing import Manager
from PyQt5.QtGui import QPixmap
from PyQt5.QtCore import Qt, QThread, QTimer, QObject, pyqtSignal
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QHBoxLayout, QVBoxLayout, QComboBox, QSlider, QLabel, QPushButton, QProgressBar
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as Figcan
import matplotlib.pyplot as plt
import backend

plt.rcParams['axes.xmargin'] = 0

class SimWorker(QObject):
    finished = pyqtSignal(object, object, object)
    def __init__(self, n_photons, batch_size, E, detectors, progress_value):
        super().__init__()
        self.n_photons = n_photons
        self.batch_size = batch_size
        self.E = E
        self.detectors = detectors
        self.progress_value = progress_value

    def run(self):
        result = backend.simulate(self.n_photons, self.batch_size, self.E, self.detectors, progress_value=self.progress_value)
        self.finished.emit(*result)

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Radiation Spectroscopy Sim")
        self.running = False
        self.patience = 1
        self.patient = None
        self.batch_prog = None
        self.figure_canvas = None
        
        #defaults
        self.E = 662 #keV
        self.batch_size = 100_000
        self.n_photons = 10_000_000

        self.manager = None
        self.progress_value = None
        self.main_layout = QHBoxLayout()

        #left controls layout
        self.left_widget = QWidget()
        self.left_layout = QVBoxLayout()
        self.left_widget.setLayout(self.left_layout)

        #First box and text
        photon_drop = QComboBox()
        photon_drop.addItems(["1_000_000", "10_000_000", "100_000_000"])
        photon_drop.setCurrentText("10_000_000")
        photon_drop.currentTextChanged.connect(self.n_photons_drop)
        photon_drop.setFixedHeight(30)
        photons_label = QLabel("Number of photons (Default = 10,000,000):")
        photons_label.setFixedHeight(30)

        #Second box and text
        batch_select = QComboBox()
        batch_select.addItems(["1_000", "10_000", "100_000"])
        batch_select.setCurrentText("100_000")
        batch_select.currentTextChanged.connect(self.batch)
        batch_select.setFixedHeight(30)
        batches_label = QLabel("Number of photons per batch (Default = 100,000):")
        batches_label.setFixedHeight(30)

        #Energy slider
        self.energy = QSlider(Qt.Horizontal)
        self.energy.setRange(1, 1000) #changed the lowest value to 1 because 0 crashes the program -K
        self.energy.setValue(662)
        self.energy.setFixedHeight(30)
        self.energy.valueChanged.connect(self.energy_update)

        #Energy label
        self.E_box = QLabel(f"Energy: {self.E} keV (Default = 662 keV)")
        self.E_box.setFixedHeight(30)

        #Start button
        start_button = QPushButton("Start")
        start_button.clicked.connect(self.main) #This runs the main that was from the backend
        start_button.setFixedHeight(30)


        #:broken_heart::wilted_rose:
        self.left_layout.addWidget(photons_label)
        self.left_layout.addWidget(photon_drop)
        self.left_layout.addWidget(batches_label)
        self.left_layout.addWidget(batch_select)
        self.left_layout.addWidget(self.E_box)
        self.left_layout.addWidget(self.energy)
        self.left_layout.addWidget(start_button)
        self.batch_prog = QProgressBar()
        self.left_layout.addWidget(self.batch_prog)
        self.left_layout.addStretch()
        #right widget placeholder for matplotlib
        self.right_widget = QWidget()
        self.right_layout = QVBoxLayout()
        self.right_widget.setLayout(self.right_layout)
        #combining into the main layout
        self.main_layout.addWidget(self.left_widget, 1)
        self.main_layout.addWidget(self.right_widget, 4)
        container = QWidget()
        container.setLayout(self.main_layout)
        self.setCentralWidget(container)
        self.progress_timer = QTimer()
        self.progress_timer.timeout.connect(self.poll_progress)

    def n_photons_drop(self, n):
        self.n_photons = int(n.replace("_", ""))

    def display_fig(self, fig):
        if self.figure_canvas:
            self.right_layout.removeWidget(self.figure_canvas)
            self.figure_canvas.deleteLater()
        self.figure_canvas = Figcan(fig)
        self.right_layout.addWidget(self.figure_canvas)

    def batch(self, b):
        self.batch_size = int(b.replace("_", ""))

    def energy_update(self, E):
        self.E = E
        self.E_box.setText(f"Energy: {E} keV (Default = 662 keV)")

    def poll_progress(self):
        if self.progress_value:
            val = self.progress_value.value
            self.batch_prog.setValue(val)
            if val >= 100:
                self.batch_prog.setValue(100)

    def main(self):
        if self.running:
            if self.patience != 0:
                self.patience = 0
                if self.patient is None:
                    self.patient = QLabel("Be Patient")
                    self.left_layout.addWidget(self.patient)
            return


        # my head hurts -K
        self.running = True
        self.patience = 1
        self.batch_prog.setValue(0)
        self.manager = Manager()
        self.progress_value = self.manager.Value('i', 0)
        self.progress_timer.start(100)
        detectors = backend.get_detectors()
        self.thread = QThread()
        self.worker = SimWorker(self.n_photons, self.batch_size, self.E,detectors, self.progress_value)
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.run)
        self.worker.finished.connect(self.simulation_done)
        self.worker.finished.connect(self.thread.quit)
        self.worker.finished.connect(self.worker.deleteLater)
        self.thread.finished.connect(self.thread.deleteLater)
        self.thread.start()

    def simulation_done(self, detected_counts, energies, total):
        self.progress_timer.stop()
        self.batch_prog.setValue(100)
        for idx, count in enumerate(detected_counts, start=1):
            print(total[idx - 1] / count * 100, "% Detector efficiency")
            print(f"detected {count} photons")

        fig = backend.plot_spectra(energies, bins=1024, energy_range=(0, 1024))#removed calling config so that it can be easier edited at this point by the user
        self.display_fig(fig)
        self.running = False

        if self.patient:
            self.left_layout.removeWidget(self.patient)
            self.patient.deleteLater()
            self.patient = None

    def closeEvent(self, event):
        raise SystemExit(0)