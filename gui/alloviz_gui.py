import sys
from PyQt5.QtWidgets import *

sys.path.append("gui")
# pyuic5 alloviz_main.ui >alloviz_main_ui.py
from alloviz_main_ui import Ui_MainWindow

from PyQt5.uic import loadUi

# See https://pypi.org/project/QtPy/
#from qtpy.QtWidgets import *
#from qtpy.uic import loadUi


class AlloVizWindow(QMainWindow, Ui_MainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setupUi(self)
        self.connectSignalsSlots()

    def connectSignalsSlots(self):
        self.actionQuit.triggered.connect(self.close)
        self.actionAbout.triggered.connect(self.about)


    def about(self):
        QMessageBox.about(
            self,
            "About Sample Editor",
            "<p>A sample text editor app built with:</p>"
            "<p>- PyQt</p>"
            "<p>- Qt Designer</p>"
            "<p>- Python</p>",
        )

if __name__ == "__main__":
    app = QApplication(sys.argv)
    win = AlloVizWindow()
    win.show()
    sys.exit(app.exec())


