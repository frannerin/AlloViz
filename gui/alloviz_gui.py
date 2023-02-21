import sys
import os
import socket
import logging
import hashlib
import json
from datetime import datetime

from PyQt5.QtWidgets import *
from PyQt5 import QtCore
from PyQt5.uic import loadUi

import pandas as pd

sys.path.append("gui")

# pyuic5 alloviz_mainwindow.ui >alloviz_mainwindow_ui.py
from alloviz_mainwindow_ui import Ui_MainWindow

import AlloViz

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(pathname)s %(levelname)s %(funcName)s %(message)s")

# See https://pypi.org/project/QtPy/
# from qtpy.QtWidgets import *
# from qtpy.uic import loadUi

# import pickle, dill
# pickle.dump(prot, open("pickle.pk","wb"))


_HOST = "localhost"
_PORT = 9990
_total_progressbar_steps = 10
_pickle_me = None

def md5sum_file(fn):
    with open(fn, "rb") as f:
        file_hash = hashlib.md5()
        while chunk := f.read(8192):
            file_hash.update(chunk)
    return file_hash.digest(), file_hash.hexdigest()

def mybreakpoint():
    from PyQt5.QtCore import pyqtRemoveInputHook, pyqtRestoreInputHook
    pyqtRemoveInputHook()
    breakpoint()
    pyqtRestoreInputHook()

def residueNumber(rn):
    """GLU:27 -> 27"""
    rnn = rn.split(":")
    return rnn[1]

""" def to_tcl_list(pylist):
    tcl_list = "{"
    for item in pylist:
        if type(item) is list:
            tcl_list += to_tcl_list(item)
        else:
            item = str(item).replace("{", "\\{").replace("}", "\\}").replace(" ", "\\ ")
            tcl_list += str(item) + " "
    tcl_list = tcl_list.strip() + "} "
    return tcl_list
 """

# TODO possibly move as inner class.
class UiStep(object):
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
            logging.error("...with exception")
            from traceback import print_tb
            print_tb(sys.exc_info()[2])
            self.obj.critical(self.msg, exc_value)
            return True


class AlloVizWindow(QMainWindow):
    updateProgress=QtCore.pyqtSignal(int)

    def __init__(self, parent=None):
        super().__init__(parent)

        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.ui.progressBar.setRange(0, _total_progressbar_steps)
        # self.ui.statusbar.addWidget(QLabel("Prova"))
        self.ui.statusbar.showMessage("Ready")
        self.fillMethodsTree()
        self.setupHistoryActions()
        self.connectSignalsSlots()

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

    def setupHistoryActions(self):
        self.ui.historyWidget.addActions([
            self.ui.actionShow_Calculation_Parameters,
            self.ui.actionViz_Analysis_Results,
            self.ui.actionOpen_Folder,
            self.ui.actionExport_Table,
        ])

    def connectSignalsSlots(self):
        self.ui.actionQuit.triggered.connect(self.close)
        self.ui.actionAbout.triggered.connect(self.showAboutDialog)
        self.ui.actionAlloViz_Homepage.triggered.connect(self.openDocumentationURL)

        self.ui.runButton.clicked.connect(self.runAnalysis)
        self.ui.methodTree.itemSelectionChanged.connect(self._updateRunButtonState)

        self.ui.actionShow_Calculation_Parameters.triggered.connect(self.historyShowCalculationParameters)
        self.ui.actionViz_Analysis_Results.triggered.connect(self.historyVizAnalysisResults)
        self.ui.actionOpen_Folder.triggered.connect(self.historyOpenFolder)
        self.ui.actionExport_Table.triggered.connect(self.historySaveAs)

        # https://stackoverflow.com/questions/50104163/update-pyqt-gui-from-a-python-thread
        self.updateProgress.connect(self.ui.progressBar.setValue)


    def historyOpenFolder(self):
        logging.info("TODO called")

    def historySaveAs(self):
        logging.info("TODO called")

    def historyShowCalculationParameters(self):
        QMessageBox.information(self,
            "Calculation Parameters",
            "TODO")

    def historyVizAnalysisResults(self):
        sl = self.ui.historyWidget.selectedItems()
        if len(sl) != 1:
            return
        uidata = sl[0].data(QtCore.Qt.UserRole)
        self.visualize(uidata)


    def _updateRunButtonState(self):
        m = self._getUiMethod()
        self.ui.runButton.setEnabled(m is not None)

    def openDocumentationURL(self):
        from PyQt5.QtCore import QUrl
        from PyQt5.QtGui import QDesktopServices
        QDesktopServices.openUrl(QUrl('https://alloviz.readthedocs.io/en/latest/?badge=latest'))

    def showAboutDialog(self):
        QMessageBox.about(
            self,
            "About the AlloViz GUI",
            "<p>The AlloViz Graphical User Interface</p>"
            "<p>A Python package to interactively compute, analyze and visualize protein allosteric communication (residue interaction) networks and delta-networks.</p>"
            "<p>Authors: Francho Nerin, Jana Selent, Toni Giorgino.</p>"
            '<p>Source code: <a href="https://github.com/frannerin/AlloViz">github.com/frannerin/AlloViz</a>.</p>'
            '<p>Please cite: <a href="https://github.com/frannerin/AlloViz">github.com/frannerin/AlloViz</a>.</p>'
        )

    def critical(self, stepname, message):
        QMessageBox.critical(
            self,
            "Error",
            f"Error while executing step `{stepname}':<br><br>{message}"
        )

    def flattenIndex(self,idx):
        """Flatten 1D or 2D indices as list"""
        idf = idx.to_frame()
        ncol = len(idf.columns)
        r = []
        for c in range(ncol):
            r.extend(idf[c].to_list())
        return r


    def checkVMDTopologyConformity(self, asel, idx):
        """Ensure that each index matches exactly one CA per result"""
        idx_flattened = self.flattenIndex(idx)
        idx_uq = list(set(idx_flattened))
        idx_uq_llist = [r.split(":") for r in idx_uq]
        r = self.doVMDcall("::alloviz::check_vmd_topology_conformity", asel, idx_uq_llist)
        return r

    def doVMDcall(self, fcn, *args):
        """Automatically serializes (via json) the arguments, then calls VMD"""
        jargs = ["::alloviz::jsonwrap", fcn]
        for a in args:
            # al = list(a)  # in case a ndarray or so
            jargs.append("{"+json.dumps(a)+"}")
        tcl = " ".join(jargs)
        r=self.sendVMDCommand(tcl)
        return r


    def sendVMDCommand(self, cmd):
        logging.info("sendVMDCommand sending: " + cmd)
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                s.connect((_HOST, _PORT))
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

    def visualize(self, uidata):
        el = uidata["elements"]
        if el=="nodes":
            self.visualizeNodes(uidata["asel"], 
                                uidata["analysis_shown"])
        elif el=="edges":
            self.visualizeEdges(uidata["asel"], 
                                uidata["analysis_shown"])
        else:
            logging.error("Should not happen")


    def visualizeNodes(self, asel, data):
        # Assumes that data is a Series
        rnl = [residueNumber(x) for x in data.index]
        rvl = data.values.tolist()
        self.doVMDcall("::alloviz::visualize_nodes", asel, rnl, rvl)


    def visualizeEdges(self, asel, data):
        # Assumes that data is a Series
        data=data.sort_values()
        data.to_csv("/tmp/debugme.csv")
        dif = data.index.to_frame()
        r1l = [residueNumber(x) for x in dif[0]]
        r2l = [residueNumber(x) for x in dif[1]]
        rvl = data.values.tolist()
        self.doVMDcall("::alloviz::visualize_edges", asel, r1l, r2l, rvl)


    def _getUiMethod(self):
        method = self.ui.methodTree.selectedItems()
        if len(method) != 1:
            return None
        kw = method[0].data(self._keyword_column, 0)
        return kw

    def _getUiFilters(self):
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
        if not flist:
            flist=["All"]   # not to be confused with "all" :(
        logging.info(f"Filters: {flist}, kwargs {fargs}") 
        return flist, fargs
    
    def _getUiAnalysisType(self):
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

    def _getUiVisualizeOptions(self):
        hide_fraction=float(self.ui.edit_hide_fraction.text())
        return hide_fraction

    def _increaseProgressBar(self):
        pbar = self.ui.progressBar
        i = pbar.value()
        logging.info(f"Updating progressbar {i} -> {i+1}")
        self.updateProgress.emit(i+1)
        #pbar.repaint()
        self.app.processEvents()

    def _showMessage(self, msg):
        logging.info(msg)
        self.ui.statusbar.showMessage(msg)
        self.app.processEvents()

    # TODO: long-running calcs in threads, interruptible
    def runAnalysis(self):
        pbar = self.ui.progressBar
        ht = self.ui.historyWidget
        asel = self.ui.atomselEdit.text()
        method = self._getUiMethod()
        filters = self._getUiFilters()

        pbar.setValue(0)
        with UiStep("Dumping trajectory", self):
            try:
                bn = self.sendVMDCommand(f"::alloviz::dump_trajectory {{{asel}}}")
                pdbfile = bn + ".pdb"
                psffile = bn + ".psf"
                dcdfile = bn + ".dcd"
                testmode = False
            except:
                logging.warning("Cannot communicate with VMD, using test data under dir 117")
                pdbfile = "117/11159_dyn_117.pdb"
                psffile = "117/11160_dyn_117.psf"
                dcdfile = "117/11156_trj_117_r.dcd"
                testmode = True

        _,tmp = md5sum_file(dcdfile)
        cache_path = f"/var/tmp/alloviz_gui_{tmp}"
        logging.info(f"Using cache path {cache_path}")

        uidata={"pdbfile": pdbfile,
                "psffile": psffile,
                "dcdfile": dcdfile,
                "asel": asel,
                "method": method
        }

        with UiStep("Loading trajectory", self):
            prot = AlloViz.Protein(pdb=pdbfile, trajs=dcdfile, path=cache_path)
            uidata["prot"] = prot

        with UiStep("Calculating", self):
            prot.calculate(method)
            calc_result = getattr(prot, method) # :( 

        with UiStep("Adding getContacts", self):
            if self.ui.checkbox_GetContacts_edges.isChecked():
                prot.calculate("GetContacts")

        with UiStep("Filtering", self):
            flist, fargs = self._getUiFilters()
            # The weird syntax requires a list of lists for sequential filtering
            prot.filter(method, filterings=[flist], **fargs)
            uidata["flist"] = flist
            uidata["fargs"] = fargs
            flist_as_string = "_".join(flist) # :( (((
            filter_result = getattr(calc_result, flist_as_string)

        with UiStep("Analyzing", self):
            el, met = self._getUiAnalysisType()
            prot.analyze(method, elements=el, metrics=met)
            uidata["elements"] = el
            uidata["metrics"] = met
            if met == "raw":
                analysis_result = filter_result._filtdata.weight
            else:
                _ = getattr(filter_result, el)
                analysis_result = getattr(_, met)
            uidata["analysis_result"] = analysis_result

        with UiStep("Checking VMD residues", self):
            if not testmode:
                if not self.checkVMDTopologyConformity(asel, analysis_result.index):
                    raise Exception("Residue numbers are not unique in the selection.")

        with UiStep("Transferring data to VMD", self):
            hide_fraction = self._getUiVisualizeOptions()
            hide_threshold = analysis_result.max() * hide_fraction
            analysis_shown = analysis_result.loc[analysis_result>hide_threshold]
            uidata["hide_fraction"] = hide_fraction
            uidata["analysis_shown"] = analysis_shown
            if not testmode:
                self.visualize(uidata)

        # Create item and add it to history 
        logging.info("Adding item to history")
        hname = datetime.now().strftime("%H:%M:%S") + " " + method
        hitem = QListWidgetItem(hname)
        hitem.setData(QtCore.Qt.UserRole, uidata)
        ht.addItem(hitem)

        self._showMessage("Ready")
        self.ui.progressBar.setValue(0)

  

        


if __name__ == "__main__":
    app = QApplication(sys.argv)
    win = AlloVizWindow()
    win.app = app
    win.show()
    sys.exit(app.exec())
