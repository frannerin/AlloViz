import sys
import os
import socket
import logging

from PyQt5.QtWidgets import *
from PyQt5 import QtCore
from PyQt5.uic import loadUi

import pandas as pd

sys.path.append("gui")

# pyuic5 alloviz_mainwindow.ui >alloviz_mainwindow_ui.py
from alloviz_mainwindow_ui import Ui_MainWindow

import AlloViz

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(pathname)s %(funcName)s %(levelname)s %(message)s")

# See https://pypi.org/project/QtPy/
# from qtpy.QtWidgets import *
# from qtpy.uic import loadUi

# TODO catch exception. possibly move as inner class.
class ComputeStep(object):
    def __init__(self, msg, obj, show_on_statusbar=True, increase_progressbar=True):
        self.msg = msg
        self.obj = obj
        self.show_on_statusbar = show_on_statusbar
        self.increase_progressbar = increase_progressbar

    def __enter__(self):
        if self.show_on_statusbar:
            self.obj.ui.statusbar.showMessage(self.msg)
        if self.increase_progressbar:
            self.obj._increaseProgressBar()
        logging.info(self.msg)
    
    def __exit__(self, exc_type, exc_value, traceback):
        logging.info("...done")
        if exc_type is not None:
            logging.info("(with exception)")
            self.obj.critical(self.msg, exc_value)
            return True


class AlloVizWindow(QMainWindow):
    updateProgress=QtCore.pyqtSignal(int)

    def __init__(self, parent=None):
        super().__init__(parent)

        self.HOST = "localhost"
        self.PORT = 9990
        self._total_progressbar_steps = 5

        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.ui.progressBar.setRange(0, self._total_progressbar_steps)
        # self.ui.statusbar.addWidget(QLabel("Prova"))
        self.ui.statusbar.showMessage("Ready")
        self.fillMethodsTree()
        self.setupHistoryWidget()

        self.connectSignalsSlots()

    def setupHistoryWidget(self):
        self.ui.historyWidget.addActions([self.ui.actionSaveAs])

    def connectSignalsSlots(self):
        self.ui.actionQuit.triggered.connect(self.close)
        self.ui.actionAbout.triggered.connect(self.about)
        self.ui.runButton.clicked.connect(self.runAnalysis)
        self.ui.methodTree.itemSelectionChanged.connect(self._updateRunButtonState)

        self.ui.actionSaveAs.triggered.connect(self.saveas)

        # https://stackoverflow.com/questions/50104163/update-pyqt-gui-from-a-python-thread
        self.updateProgress.connect(self.ui.progressBar.setValue)


    def saveas(self):
        logging.info("SAVEAS called")


    def _updateRunButtonState(self):
        m = self._getMethod()
        self.ui.runButton.setEnabled(m is not None)


    def fillMethodsTree(self):
        wlist = AlloViz.AlloViz.info.wrappers
        df = pd.DataFrame(wlist.values())
        df[4] = wlist.keys()
        df.columns = ["Quantity", "Software", "Correlation metric", "Object", "Keyword"]

        # df.set_index(["Metric","Quantity","Object"])
        #  Rearrange the columns according to the desired nesting
        df = df[["Quantity", "Object", "Correlation metric", "Software", "Keyword"]]
        self._keyword_column = 4

        # Remove gRINN because it requires ff parameters
        df=df.loc[df.Software != "gRINN"]

        # tree = self.findChild(QTreeWidget,"methodTree")
        tree = self.ui.methodTree
        tree.setColumnCount(len(df.columns))
        tree.setHeaderLabels(df.columns)


        # Nest grouping the first 2 levels by unique value.
        items = []
        ditto = "..."  # "〃"
        for lev1 in df.iloc[:, 0].unique():
            lev1_item = QTreeWidgetItem([lev1])
            lev1_item.setFlags(lev1_item.flags() & ~QtCore.Qt.ItemIsSelectable)
            # lev1_item.setFirstColumnSpanned(True)  ## would be nice but does nothing
            items.append(lev1_item)
            df1 = df.loc[df.iloc[:, 0] == lev1]
            for lev2 in df1.iloc[:, 1].unique():
                lev2_item = QTreeWidgetItem([ditto, lev2])
                lev2_item.setFlags(lev2_item.flags() & ~QtCore.Qt.ItemIsSelectable)
                lev1_item.addChild(lev2_item)
                df2 = df1.loc[df.iloc[:, 1] == lev2]
                for i, leaf in df2.iterrows():
                    leaf_list = list(leaf)
                    leaf_list[0] = ditto
                    leaf_list[1] = ditto
                    leaf_item = QTreeWidgetItem(leaf_list)
                    lev2_item.addChild(leaf_item)
        tree.insertTopLevelItems(0, items)

        # tree.setColumnWidth(0,200)
        tree.resizeColumnToContents(0)

        # for i in range(len(df.columns)):
        #    tree.resizeColumnToContents(i)

    def about(self):
        QMessageBox.about(
            self,
            "About the AlloViz GUI",
            "<p>The AlloViz Graphical User Interface</p>"
            "<p>Authors: Francho Nerin, Jana Selent, Toni Giorgino</p>"
            '<p>Source code: <a href="https://github.com/frannerin/AlloViz">github.com/frannerin/AlloViz</a></p>',
        )

    def critical(self, stepname, message):
        QMessageBox.critical(
            self,
            "Error",
            f"Error while executing step `{stepname}':<br><br>{message}"
        )

    def sendVMDCommand(self, cmd):
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                s.connect((self.HOST, self.PORT))
                s.sendall(str.encode(cmd + "\n"))
                data = s.recv(4096)
        except Exception as e:
            logging.warning(
                f"Could not connect to VMD ({self.HOST}:{self.PORT}): {e}.\nMake sure the client component is running."
            )
            raise e

        ret = data.decode().strip()
        logging.info("sendVMDCommand got " + ret)
        return ret

    def _getMethod(self):
        method = self.ui.methodTree.selectedItems()
        if len(method) != 1:
            return None
        kw = method[0].data(self._keyword_column, 0)
        return kw

    def _getFilters(self):
        flist = []
        fargs = {}
        if self.ui.checkbox_GetContacts_edges.isChecked():
            flist.append("GetContacts_edges")
            fargs["GetContacts_threshold"]=float(self.ui.edit_GetContacts_threshold.text())
        if self.ui.checkbox_No_Sequence_Neighbors.isChecked():
            flist.append("No_Sequence_Neighbors")
            fargs["Sequence_Neighbor_distance"]=int(self.ui.edit_Sequence_Neighbor_distance.text())
        if self.ui.checkbox_Spatially_distant.isChecked():
            flist.append("Spatially_distant")
            fargs["Interresidue_distance"]=float(self.ui.edit_Interresidue_distance.text())
        if self.ui.checkbox_GPCR_Interhelix.isChecked():
            flist.append("GPCR_Interhelix")
        logging.info(f"Filters: {flist}, kwargs {fargs}") 
        return flist, fargs
    
    def _getAnalysis(self):
        if self.ui.anEdgeBtwCheck.isChecked():
            el, met = "edges", "btw"
        elif self.ui.anEdgeCurrentCheck.isChecked():
            el, met = "edges", "cfb"
        elif self.ui.anEdgeRawCheck.isChecked():
            el, met = "edges", "raw" # ?????
        elif self.ui.anNodeBtwCheck.isChecked():
            el, met = "nodes", "btw"
        elif self.ui.anNodeCurrentCheck.isChecked():
            el, met = "nodes", "cfb"
        else:
            logging.error("No radio button selected, should not happen")
        logging.info(f"Analysis: {el}, {met}")
        return el, met

    # Should work but doesn't update. Not with setValue, neither with signals
    def _increaseProgressBar(self):
        pbar = self.ui.progressBar
        i = pbar.value()
        logging.info(f"Current progressbar is {i}")
        self.updateProgress.emit(i+1)
        #pbar.repaint()
        self.app.processEvents()

    def _showMessage(self, msg):
        logging.info(msg)
        self.ui.statusbar.showMessage(msg)
        self.app.processEvents()

    def runAnalysis(self):
        pbar = self.ui.progressBar
        ht = self.ui.historyWidget
        asel = self.ui.atomselEdit.text()
        method = self._getMethod()
        filters = self._getFilters()

        logging.info(f"FOCUS ht: {ht.hasFocus()}")
        logging.info(f"FOCUS m: {self.ui.methodTree.hasFocus()}")

        logging.info(f"Run clicked: {asel}, {method}")

        pbar.setValue(0)
        with ComputeStep("Dumping trajectory", self):
            try:
                bn = self.sendVMDCommand(f"::alloviz::dump_trajectory {{{asel}}}")
                pdbfile = bn + ".pdb"
                dcdfile = bn + ".dcd"
            except:
                logging.warning("Cannot communicate with VMD, using test data under dir 117")
                pdbfile = "../117/11159_dyn_117.pdb"
                dcdfile = "../117/11157_trj_117.xtc"

        with ComputeStep("Loading trajectory", self):
            prot = AlloViz.Protein(pdb=pdbfile, trajs=dcdfile)

        self._increaseProgressBar()
        with ComputeStep("Calculating", self):
            prot.calculate(method)

        self._increaseProgressBar()
        with ComputeStep("Adding getContacts", self):
            if self.ui.checkbox_GetContacts_edges.isChecked():
                prot.calculate("GetContacts")

        self._increaseProgressBar()
        with ComputeStep("Filtering", self):
            flist, fargs = self._getFilters()
            # The weird syntax requires a list of lists for sequential filtering
            prot.filter("all", filterings=[flist], **fargs)


        self._increaseProgressBar()
        with ComputeStep("Analyzing", self):
            el, met = self._getAnalysis()
            # The weird syntax requires a list of lists for sequential filtering
            prot.analyze("all", elements=el, metrics=met)

        self._showMessage("Ready")

        # Add item to history 
        hitem = QListWidgetItem(method)
        hitem.setData(QtCore.Qt.UserRole, prot)
        ht.addItem(hitem)

        


if __name__ == "__main__":
    app = QApplication(sys.argv)
    win = AlloVizWindow()
    win.app = app
    win.show()
    sys.exit(app.exec())
