#GUI
import sys
from PyQt5 import QtCore as Qtcre, QtWidgets as Qtwdgt, QtGui as Qtgui
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as Figcan

class MyWidget(Qtwdgt.QWidget):
    def __init__(self, fig):
        super().__init__()
        self.text = Qtwdgt.QLabel("hello world", alignment=Qtcre.Qt.AlignCenter)
        self.toolbar = Figcan(fig)
        
        self.layout = Qtwdgt.QVBoxLayout(self)
        self.layout.addWidget(self.text)
        self.layout.addWidget(self.toolbar)

