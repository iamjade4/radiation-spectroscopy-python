from frontend import MainWindow
from PyQt5.QtWidgets import QApplication
import sys

app = QApplication(sys.argv)
window = MainWindow()
window.resize(1000, 800)
window.show()
app.exec()
