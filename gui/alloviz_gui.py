import sys
import os
from tkinter.tix import COLUMN
from PyQt5.QtWidgets import *
import pandas as pd

# os.environ["DYLD_LIBRARY_PATH"]="/opt/homebrew/opt/sqlite/lib/"

sys.path.append("gui")
# pyuic5 alloviz_main.ui >alloviz_main_ui.py
from alloviz_main_ui import Ui_MainWindow

import AlloViz

from PyQt5.uic import loadUi

# See https://pypi.org/project/QtPy/
#from qtpy.QtWidgets import *
#from qtpy.uic import loadUi


class AlloVizWindow(QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.fillMethodsTree()
        self.connectSignalsSlots()

    def connectSignalsSlots(self):
        self.ui.actionQuit.triggered.connect(self.close)
        self.ui.actionAbout.triggered.connect(self.about)

    def fillMethodsTree(self):
        wlist=AlloViz.AlloViz.info.wrappers
        df = pd.DataFrame(wlist.values())
        df[4]=wlist.keys()
        df.columns = ["Quantity", "Software", "Metric", "Object", "Code"]

        tree = self.findChild(QTreeWidget,"methodTree")
        tree.setColumnCount(5)
        tree.setHeaderLabels(["Quantity", "Metric", "Software", "Object", "Code"])

        items = []
        for qty in df.Quantity.unique():
            qty_item = QTreeWidgetItem([qty])
            items.append(qty_item)
            dfs = df.loc[df.Quantity == qty]            
            for metr in dfs.Metric.unique():
                metr_item = QTreeWidgetItem(["", metr])
                qty_item.addChild(metr_item)
                dfsm = dfs.loc[df.Metric==metr]
                for i,poc in dfsm.iterrows():
                    poc_item = QTreeWidgetItem(["", "", poc.Package, poc.Object, poc.Code])
                    metr_item.addChild(poc_item)
        tree.insertTopLevelItems(0, items)

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


