import sys
from multiprocessing import Manager
from PyQt5.QtGui import QPixmap
from PyQt5.QtCore import Qt, QThread, QTimer, QObject, pyqtSignal
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QHBoxLayout, QVBoxLayout, QComboBox, QSlider, QLabel, QPushButton, QProgressBar, QCheckBox, QLineEdit
from PyQt5.QtGui import QIntValidator #THE VALIDATORRR
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as Figcan
from detectors.naitl import NaITl
from detectors.si import Si
import time
import matplotlib.pyplot as plt
import numpy as np
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

class Detector(QWidget): #Window to add new detectors
    def __init__(self, parent):
        super(Detector, self).__init__(parent=None)
        self.parent = parent
        self.cyl_defaults = ['NaITl', 600, 100, 0, 6, 1, 0, 60, 120]
        self.sq_defaults = ['Si', 120, 120, 1000, 1, 0, 0, 600, 0, 0]
        #self.detector = NaITl([600, 100, 0], [self.axis_x, self.axis_y, self.axis_z]/(np.sqrt(self.axis_x**2 + self.axis_y**2 + self.axis_z**2)), self.radius, self.height)
        self.setWindowTitle("New Detector")
        self.cyl_layout = QVBoxLayout() #Layout for when a cylindrical detector is being created
        self.cyl_widget = QWidget()
        self.cyl_widget.setLayout(self.cyl_layout)
        self.square_layout = QVBoxLayout()
        self.square_widget = QWidget()
        self.square_widget.setLayout(self.square_layout)
        self.layout = QVBoxLayout() #Mainlayout that the specific layouts get added to
        
        self.add_label = QLabel("Add a new detector")
        self.add_label.setFixedHeight(30)
        
        self.options = QComboBox()
        self.options.addItem('')
        self.options.addItem("NaITl")
        self.options.addItem("Si")
        self.options.setFixedHeight(30)
        self.options.setCurrentIndex(0)
        self.options.currentIndexChanged.connect(self.detector_select)
        
        #widgets for cyl_layouts
        self.radius_label = QLabel("Radius (mm): (Default = 60)")
        self.radius_label.setFixedHeight(30)
        
        self.radius_input = QLineEdit(parent=self)
        self.radius_input.setValidator(QIntValidator(bottom=1))
        self.radius_input.setFixedHeight(30)
        
        self.cyl_height_label = QLabel("Height (mm): (Default = 120)")
        self.cyl_height_label.setFixedHeight(30)
        
        self.cyl_height_input = QLineEdit(parent=self)
        self.cyl_height_input.setValidator(QIntValidator(bottom=1))
        self.cyl_height_input.setFixedHeight(30)
        
        self.axis_label = QLabel("Cylindrical axis: (Default = [6, 1, 0]). Input in X, Y, Z order.")
        self.axis_label.setFixedHeight(30)
        
        self.axis_input_x = QLineEdit(parent=self)
        self.axis_input_x.setValidator(QIntValidator())
        self.axis_input_x.setFixedHeight(30)
        self.axis_input_y = QLineEdit(parent=self)
        self.axis_input_y.setValidator(QIntValidator())
        self.axis_input_y.setFixedHeight(30)
        self.axis_input_z = QLineEdit(parent=self)
        self.axis_input_z.setValidator(QIntValidator())
        self.axis_input_z.setFixedHeight(30)
        
        self.base_label = QLabel("Position of the detector: (Default = [600, 100, 0]. Input in X, Y, Z order.")
        self.base_label.setFixedHeight(30)
                                       
        self.base_input_x = QLineEdit(parent=self)
        self.base_input_x.setValidator(QIntValidator())
        self.base_input_x.setFixedHeight(30)
        self.base_input_y = QLineEdit(parent=self)
        self.base_input_y.setValidator(QIntValidator())
        self.base_input_y.setFixedHeight(30)
        self.base_input_z = QLineEdit(parent=self)
        self.base_input_z.setValidator(QIntValidator())
        self.base_input_z.setFixedHeight(30)

        #widgets for square layouts
        self.sq_height_label = QLabel("Height (mm): (Default = 120)")
        self.sq_height_label.setFixedHeight(30)
        
        self.sq_height_input = QLineEdit(parent=self)
        self.sq_height_input.setValidator(QIntValidator(bottom=1))
        self.sq_height_input.setFixedHeight(30)
        
        self.width_label = QLabel("Width (mm): (Default = 120)")
        self.width_label.setFixedHeight(30)
        
        self.width_input = QLineEdit(parent=self)
        self.width_input.setValidator(QIntValidator(bottom=1))
        self.width_input.setFixedHeight(30)
        
        self.depth_label = QLabel("Thickness (um): (Default = 1000)")
        self.depth_label.setFixedHeight(30)
        
        self.depth_input = QLineEdit(parent=self)
        self.depth_input.setValidator(QIntValidator(bottom=1))
        self.depth_input.setFixedHeight(30)
        
        self.normal_label = QLabel("Direction the detector is facing: (Default = [1, 0, 0]). Input in X, Y, Z order.")
        self.normal_label.setFixedHeight(30)
        
        self.normal_input_x = QLineEdit(parent=self)
        self.normal_input_x.setValidator(QIntValidator())
        self.normal_input_x.setFixedHeight(30)
        self.normal_input_y = QLineEdit(parent=self)
        self.normal_input_y.setValidator(QIntValidator())
        self.normal_input_y.setFixedHeight(30)
        self.normal_input_z = QLineEdit(parent=self)
        self.normal_input_z.setValidator(QIntValidator())
        self.normal_input_z.setFixedHeight(30)
        
        self.pos_label = QLabel("Position of the centre front of the detector (Default = [600, 0, 0]). Input in X, Y, Z order.")
        self.pos_label.setFixedHeight(30)
        
        self.pos_input_x = QLineEdit(parent=self)
        self.pos_input_x.setValidator(QIntValidator())
        self.pos_input_x.setFixedHeight(30)
        self.pos_input_y = QLineEdit(parent=self)
        self.pos_input_y.setValidator(QIntValidator())
        self.pos_input_y.setFixedHeight(30)
        self.pos_input_z = QLineEdit(parent=self)
        self.pos_input_z.setValidator(QIntValidator())
        self.pos_input_z.setFixedHeight(30)
        
        self.add = QPushButton("Add Detector")
        self.add.setFixedHeight(30)
        self.add.clicked.connect(self.add_detector)

        self.layout.addWidget(self.add_label)
        self.layout.addWidget(self.options) #main layout regardless of detector
        
        
        self.cyl_layout.addWidget(self.radius_label)
        self.cyl_layout.addWidget(self.radius_input)
        self.cyl_layout.addWidget(self.cyl_height_label)
        self.cyl_layout.addWidget(self.cyl_height_input)
        self.cyl_layout.addWidget(self.axis_label)
        self.cyl_layout.addWidget(self.axis_input_x)
        self.cyl_layout.addWidget(self.axis_input_y)
        self.cyl_layout.addWidget(self.axis_input_z)
        self.cyl_layout.addWidget(self.base_label)
        self.cyl_layout.addWidget(self.base_input_x)
        self.cyl_layout.addWidget(self.base_input_y)
        self.cyl_layout.addWidget(self.base_input_z)
        
        self.square_layout.addWidget(self.sq_height_label)
        self.square_layout.addWidget(self.sq_height_input)
        self.square_layout.addWidget(self.width_label)
        self.square_layout.addWidget(self.width_input)
        self.square_layout.addWidget(self.depth_label)
        self.square_layout.addWidget(self.depth_input)
        self.square_layout.addWidget(self.normal_label)
        self.square_layout.addWidget(self.normal_input_x)
        self.square_layout.addWidget(self.normal_input_y)
        self.square_layout.addWidget(self.normal_input_z)
        self.square_layout.addWidget(self.pos_label)
        self.square_layout.addWidget(self.pos_input_x)
        self.square_layout.addWidget(self.pos_input_y)
        self.square_layout.addWidget(self.pos_input_z)
        
        self.layout.addWidget(self.cyl_widget)
        self.cyl_widget.setVisible(0)
        self.layout.addWidget(self.square_widget)
        self.square_widget.setVisible(0)
        self.layout.addWidget(self.add)
        self.parameters = self.cyl_defaults
        self.setLayout(self.layout)


    def detector_select(self):
        index = self.options.currentIndex()
        if index == 1:
            self.square_widget.setVisible(0)
            self.cyl_widget.setVisible(1)
        elif index == 2:
            self.cyl_widget.setVisible(0)
            self.square_widget.setVisible(1)

    def add_detector(self):
        if self.options.currentIndex() == 1:
            self.inputs = [self.options.currentText(), self.base_input_x.text(), self.base_input_y.text(), self.base_input_z.text(), self.axis_input_x.text(), self.axis_input_y.text(), self.axis_input_z.text(), self.radius_input.text(), self.cyl_height_input.text()]
            self.parameters = self.cyl_defaults
        elif self.options.currentIndex() ==2:
            print("yerp")
            self.inputs = [self.options.currentText(), self.sq_height_input.text(), self.width_input.text(), self.depth_input.text(), self.normal_input_x.text(), self.normal_input_y.text(), self.normal_input_z.text(), self.pos_input_x.text(), self.pos_input_y.text(), self.pos_input_z.text()]
            self.parameters = self.sq_defaults
        for i in range(len(self.inputs)):
            if self.inputs[i] != '':
                self.parameters[i] = self.inputs[i]
            #VV Just putting NaITl manually here for now
        if self.parameters[0] == 'NaITl':
            self.detector = NaITl([int(self.parameters[1]), int(self.parameters[2]), int(self.parameters[3])], [int(self.parameters[4]), int(self.parameters[5]), int(self.parameters[6])]/(np.sqrt(int(self.parameters[4])**2 + int(self.parameters[5])**2 + int(self.parameters[6])**2)), int(self.parameters[7]), int(self.parameters[8]))
        elif self.parameters[0] == 'Si':
            self.detector = Si(int(self.parameters[1]), int(self.parameters[2]), int(self.parameters[3]), [int(self.parameters[4]), int(self.parameters[5]), int(self.parameters[6])]/(np.sqrt(int(self.parameters[4])**2 + int(self.parameters[5])**2 + int(self.parameters[6])**2)), [int(self.parameters[7]), int(self.parameters[8]), int(self.parameters[9])])
        self.parent.add_detector_func()
        self.close()
class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)
        self.setWindowTitle("Radiation Spectroscopy Sim")
        self.running = False
        self.patience = 1
        self.patient = None
        self.batch_prog = None
        self.figure_canvas = None
        self.detect_num = 0
        
        #defaults
        self.E = 662 #keV
        self.batch_size = 100_000
        self.n_photons = 10_000_000
        self.detectors = [
            NaITl([600, 100, 0], [6, 1, 0]/(np.sqrt(37)), 60, 120),  #cylinder with a radius of 6cm, height of 12cm, axis facing the origin
            Si(120, 120, 1000, [1, 0, 0], [600, 0, 000])       # 12x12x0.1cm cuboidal detector
        ]

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

        #Energy label
        self.E_box = QLabel(f"Energy: {self.E} keV (Default = 662 keV)")
        self.E_box.setFixedHeight(30)

        #Energy slider
        self.energy_slider = QSlider(Qt.Horizontal)
        self.energy_slider.setRange(1, 1000)
        self.energy_slider.setValue(self.E)
        self.energy_slider.setFixedHeight(30)
        self.energy_slider.valueChanged.connect(self.energy_update)

        #Self-validating integer input box
        self.energy_input = QLineEdit()
        self.energy_input.setValidator(QIntValidator(1, 1000)) #so the interesting chungus quirk about this is that if it isnt between 1-1000 it just uses the previous number
        self.energy_input.setFixedHeight(30)
        self.energy_input.setVisible(False)
        self.energy_input.setPlaceholderText("Energy in keV (1–1000)")
        self.energy_input.returnPressed.connect(self.energy_input_update)

        #Checkbox for input toggling
        self.use_text_input = QCheckBox("Use input box instead")
        self.use_text_input.setFixedHeight(30)
        self.use_text_input.stateChanged.connect(self.toggle_energy_input)
        
        self.add_detector_label = QLabel("""If no detectors are added, default detectors are:
                                         NaITl([600, 100, 0], [6, 1, 0]/(np.sqrt(37)), 60, 120)
                                         Si(100, 0, 100, 60, 60, 2)""")
        self.add_detector_label.setFixedHeight(60)
        
        self.add_detector = QPushButton("Add detector")
        self.add_detector.clicked.connect(self.detector)
        self.add_detector.setFixedHeight(30)
        
        self.detector_lim = QLabel("You have reached the detector limit of 4.")
        self.detector_lim.setFixedHeight(30)
        self.detector_lim.hide()

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
        self.left_layout.addWidget(self.energy_slider)
        self.left_layout.addWidget(self.energy_input)
        self.left_layout.addWidget(self.use_text_input)
        self.left_layout.addWidget(self.add_detector_label)
        self.left_layout.addWidget(self.add_detector)
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
        if self.use_text_input.isChecked():
            self.energy_input.setText(str(E))

    def poll_progress(self):
        if self.progress_value:
            val = self.progress_value.value
            self.batch_prog.setValue(val)
            if val >= 100:
                self.batch_prog.setValue(100)

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

    def toggle_energy_input(self, state):
        if state == Qt.Checked:
            self.energy_slider.setVisible(False)
            self.energy_input.setVisible(True)
            self.energy_input.setText(str(self.energy_slider.value()))
        else:
            text = self.energy_input.text()
            if text.isdigit():
                self.energy_update(int(text))
            self.energy_input.setVisible(False)
            self.energy_slider.setVisible(True)

    def energy_input_update(self):
        text = self.energy_input.text()
        if text.isdigit():
            val = int(text)
            if 1 <= val <= 1000:
                self.E = val
                self.E_box.setText(f"Energy: {val} keV (Default = 662 keV)")

    def closeEvent(self, event):
        raise SystemExit(0)
    
    def detector(self):
        if self.detect_num == 4:
            self.detector_lim.show()
        else:
            self.window = Detector(self)
            self.window.show()
        
    def add_detector_func(self):
        self.detect_num += 1
        if self.detect_num == 1:
            self.detectors = [] #deletes the defaults only if it's the first time the function is called
        self.detectors.append(self.window.detector)

    #main is down here now because its faster to get to it if its at either end of the file 
        #oh ok
    def main(self):
        if self.use_text_input.isChecked():
            text = self.energy_input.text()
            if text.isdigit():
                val = int(text)
                if 1 <= val <= 1000:
                    self.energy_update(val)

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
        detectors = self.detectors
        self.thread = QThread()
        self.worker = SimWorker(self.n_photons, self.batch_size, self.E,detectors, self.progress_value)
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.run)
        self.worker.finished.connect(self.simulation_done)
        self.worker.finished.connect(self.thread.quit)
        self.worker.finished.connect(self.worker.deleteLater)
        self.thread.finished.connect(self.thread.deleteLater)
        self.thread.start()
