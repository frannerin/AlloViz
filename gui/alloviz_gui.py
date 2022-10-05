import sys
import os
from PyQt5.QtWidgets import *
import pandas as pd

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
        df.columns = ["Quantity", "Software", "Correlation metric", "Object", "Keyword"]

        #df.set_index(["Metric","Quantity","Object"])
        # Rearrange the columns according to the desired nesting
        df = df[["Correlation metric", "Quantity", "Object", "Software", "Keyword"]]

        tree = self.findChild(QTreeWidget,"methodTree")
        tree.setColumnCount(len(df.columns))
        tree.setHeaderLabels(df.columns)

        items = []
        ditto = "〃"
        for lev1 in df.iloc[:,0].unique():
            lev1_item = QTreeWidgetItem([lev1])
            items.append(lev1_item)
            df1 = df.loc[df.iloc[:,0] == lev1]            
            for lev2 in df1.iloc[:,1].unique():
                lev2_item = QTreeWidgetItem([ditto, lev2])
                lev1_item.addChild(lev2_item)
                df2 = df1.loc[df.iloc[:,1]==lev2]
                for i,leaf in df2.iterrows():
                    leaf_list = list(leaf)
                    leaf_list[0] = ditto
                    leaf_list[1] = ditto
                    leaf_item = QTreeWidgetItem(leaf_list)
                    lev2_item.addChild(leaf_item)
        tree.insertTopLevelItems(0, items)
        
        #tree.setColumnWidth(0,200)
        tree.resizeColumnToContents(0)

        #for i in range(len(df.columns)):
        #    tree.resizeColumnToContents(i)

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


