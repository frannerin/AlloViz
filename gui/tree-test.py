import sys

# https://www.pythonguis.com/faq/pyqt5-vs-pyside2/
from qtpy.QtWidgets import QApplication, QTreeWidget, QTreeWidgetItem

# https://doc.qt.io/qtforpython/tutorials/basictutorial/treewidget.html

data = {"Project A": ["file_a.py", "file_a.txt", "something.xls"],
        "Project B": ["file_b.csv", "photo.jpg"],
        "Project C": []}

app = QApplication(sys.argv)

tree = QTreeWidget()
tree.setColumnCount(2)
tree.setHeaderLabels(["Name", "Type"])

items = []
for key, values in data.items():
    item = QTreeWidgetItem([key])
    for value in values:
        ext = value.split(".")[-1].upper()
        child = QTreeWidgetItem([value, ext])
        item.addChild(child)
    items.append(item)

tree.insertTopLevelItems(0, items)

tree.show()
app.exec()
