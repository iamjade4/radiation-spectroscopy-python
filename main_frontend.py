#GUI
import sys
from PyQt5.QtCore import QSize, Qt
from PyQt5.QtWidgets import QApplication, QWidget, QMainWindow, QPushButton, QComboBox, QGridLayout, QSlider, QLabel
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as Figcan
from main_backend import main

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
        self.energy.setRange(0,1000)
        self.energy.setSingleStep(1)
        self.energy.setFixedHeight(30)
        self.energy.setValue(662)
        self.E = self.energy.value()
        self.energy.valueChanged.connect(self.energy_update)
        layout.addWidget(self.energy,5,0)
        
        #Energy label
        self.E_box = QLabel("Energy: "+str(self.E)+"keV")
        self.E_box.setFixedHeight(30)
        layout.addWidget(self.E_box,4,0)
        
        #Start button
        start_button = QPushButton("Start")
        start_button.clicked.connect(self.run_main)
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
        print(self.n_photons)

    
    def display_fig(self, fig):
        figure = Figcan(fig)
        height = self.layout.rowCount()
        self.layout.addWidget(figure,0,1,height,3)
        
    def batch(self, b):
        self.batch_size = int(b)
        print(self.batch_size)
    
    def energy_update(self, E):
        self.E_box.setText("Energy: "+str(E)+"keV")
        self.E = int(E)
        
    def run_main(self):
        fig = main(self.n_photons, self.batch_size, self.E)
        self.display_fig(fig)

app = QApplication(sys.argv)
window = MainWindow()
window.resize(1000,800)
window.show()
app.exec()  
