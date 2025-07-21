
import sys

from PyQt5.QtCore import Qt

from PyQt5.QtWidgets import (
    QGroupBox,
    QApplication,
    QCheckBox,
    QComboBox,
    QDateEdit,
    QDateTimeEdit,
    QDial,
    QDoubleSpinBox,
    QFontComboBox,
    QLabel,
    QLCDNumber,
    QLineEdit,
    QMainWindow,
    QProgressBar,
    QPushButton,
    QRadioButton,
    QSlider,
    QSpinBox,
    QTimeEdit,
    QHBoxLayout,
    QVBoxLayout,
    QFrame,
    QWidget,
)

import pandas as pd
import mrich
from mrich import print
from pathlib import Path

DIRECTORY = Path(".").resolve()
CONFIG_FILE = Path("config.yaml").resolve()
ASSEMBLIES_FILE = Path("assemblies.yaml").resolve()

class XChemAlign(QWidget):

    def __init__(self):
        super().__init__()
        self.setup_window()
        self.setup_layout()
        self.ui_config = ConfigUI(self, self.layout_config)

    def setup_window(self, width=400, height=200):
        self.setWindowTitle("XChemAlign GUI")
        self.resize(width,height)

    def setup_layout(self):
        
        self.layout = QHBoxLayout()

        
        self.layout_left = QVBoxLayout()
        self.layout_config = QVBoxLayout()
        label = QLabel("XChemAlign")
        font = label.font()
        font.setPointSize(16)   # bigger size
        font.setBold(True)      # bold text
        label.setFont(font)
        label.setAlignment(Qt.AlignCenter)
        self.layout_config.addWidget(label)
        
        label = QLabel(f"{DIRECTORY}")
        font = label.font()
        font.setFamily("Courier New")  # or "Consolas", "Monospace", etc.
        label.setFont(font)
        label.setStyleSheet("background-color:black; color:white; padding:8px")
        self.layout_config.addWidget(label)
        # self.layout_assemblies = QVBoxLayout()
        # self.layout_console = QVBoxLayout()

        self.layout.addLayout(self.layout_left)
        # self.layout.addLayout(self.layout_console)
        
        self.layout_left.addLayout(self.layout_config)
        # self.layout_left.addLayout(self.layout_assemblies)

        label = QLabel("CONFIG")
        font = label.font()
        font.setPointSize(14)
        font.setBold(True)
        label.setAlignment(Qt.AlignCenter)
        label.setFont(font)
        self.layout_config.addWidget(label)
        # self.layout_assemblies.addWidget(QLabel("ASSEMBLIES"))
        # self.layout_console.addWidget(QLabel("CONSOLE"))

        self.setLayout(self.layout)

    def on_click(self):
        self.button.setText("You clicked me!")

class ChildUI(QWidget):
    def __init__(self, parent, layout, index=None):
        super().__init__()
        self.parent = parent
        self.layout = layout
        self.index = index
        self.setup()
    
class ConfigUI(ChildUI):
    def setup(self):
        
        # target name
        self.label_target_name = QLabel("Target Name:")
        self.widget_target_name = QLineEdit()
        l = QHBoxLayout()
        l.addWidget(self.label_target_name)
        l.addWidget(self.widget_target_name)
        self.layout.addLayout(l)

        # target name
        self.label_base_dir = QLabel("Base Directory:")
        self.widget_base_dir = QLineEdit("/")
        l = QHBoxLayout()
        l.addWidget(self.label_base_dir)
        l.addWidget(self.widget_base_dir)
        self.layout.addLayout(l)

        # # inputs
        self.layout_inputs = QVBoxLayout()
        self.layout.addLayout(self.layout_inputs)
        self.ui_inputs = InputsUI(self, self.layout_inputs)

        # panddas missing ok

        # ref datasets

    def load_yaml(self, file):

        import yaml
        with open(file, "r") as stream:
            config = yaml.safe_load(stream)

        self.widget_target_name.setText(config["target_name"])
        self.widget_base_dir.setText(config["base_dir"])

        self.ui_inputs.clear()

        for d in config["inputs"]:
            self.ui_inputs.add_input(
                dir=d["dir"],
                type=d["type"],
                code_prefix=d["code_prefix"],
                code_prefix_tooltip=d["code_prefix_tooltip"],
                soakdb=d.get("soakdb",""),
            )

class InputsUI(ChildUI):

    def clear(self):
        self.inputs = pd.DataFrame(columns=["type", "code_prefix", "code_prefix_tooltip", "soakdb"])

    def setup(self):

        self.clear()

        self.label_title = QLabel("Inputs")
        font = self.label_title.font()
        font.setBold(True)
        self.label_title.setFont(font)
        self.layout.addWidget(self.label_title)
        label = QLabel("Define your inputs here:")
        self.layout.addWidget(label)

        self.button_add = QPushButton("Add input")
        self.button_add.clicked.connect(self.add_input)
        self.layout.addWidget(self.button_add)

        self.layout_inputs = QVBoxLayout()
        self.layout.addLayout(self.layout_inputs)

        self.layout_buttons = QHBoxLayout()
        
        b1 = QPushButton("SAVE")
        b1.setStyleSheet("color:green")
        def b1_action():
            raise NotImplementedError
        b1.clicked.connect(b1_action)
        self.layout_buttons.addWidget(b1)

        b2 = QPushButton("REVERT")
        b2.setStyleSheet("color:red")
        def b2_action():
            self.parent.load_yaml(CONFIG_FILE)
        b2.clicked.connect(b2_action)
        self.layout_buttons.addWidget(b2)

        
        self.layout.addLayout(self.layout_buttons)

        ###

        self.update_widgets()

    def add_input(self,
        *args,
        dir="",
        type="model_building",
        code_prefix="",
        code_prefix_tooltip="",
        soakdb="",
    ):

        df = pd.DataFrame([dict(
            dir=dir,
            type=type,
            code_prefix=code_prefix,
            code_prefix_tooltip=code_prefix_tooltip,
            soakdb=soakdb,
        )])

        self.inputs = pd.concat([self.inputs, df], ignore_index=True)

        print(self.inputs)

        self.update_widgets()

    def update_widgets(self):

        clear_layout(self.layout_inputs)

        for i,row in self.inputs.iterrows():

            groupbox = QGroupBox()

            layout = QVBoxLayout()

            groupbox.setLayout(layout)
            
            # title
            l = QHBoxLayout()
            label = QLabel(f"INPUT {i+1}")
            font = label.font()
            # font.setPointSize(14)
            font.setBold(True)
            # label.setAlignment(Qt.AlignCenter)
            label.setFont(font)
            l.addWidget(label)

            def remove_input():
                self.inputs = self.inputs.drop(i)
                self.update_widgets()

            # remove
            button = QPushButton("remove")
            button.clicked.connect(remove_input)
            button.setStyleSheet("color: red;")
            l.addWidget(button)

            layout.addLayout(l)

            # dir
            label = QLabel("dir")
            widget = QLineEdit(row["dir"])
            l = QHBoxLayout()
            l.addWidget(label)
            l.addWidget(widget)
            layout.addLayout(l)

            # type
            label = QLabel("type")
            widget = QComboBox()
            widget.addItems(["model_building", "manual"])
            index = widget.findText(row["type"])
            if index != -1:
                widget.setCurrentIndex(index)
            l = QHBoxLayout()
            l.addWidget(label)
            l.addWidget(widget)
            layout.addLayout(l)

            # code_prefix
            label = QLabel("code_prefix")
            widget = QLineEdit(row["code_prefix"])
            l = QHBoxLayout()
            l.addWidget(label)
            l.addWidget(widget)
            layout.addLayout(l)

            # code_prefix
            label = QLabel("code_prefix_tooltip")
            widget = QLineEdit(row["code_prefix_tooltip"])
            l = QHBoxLayout()
            l.addWidget(label)
            l.addWidget(widget)
            layout.addLayout(l)

            self.layout_inputs.addWidget(groupbox)

def clear_layout(layout):
    if layout is not None:
        while layout.count():
            item = layout.takeAt(0)

            widget = item.widget()
            if widget is not None:
                widget.setParent(None)
                widget.deleteLater()

            child_layout = item.layout()
            if child_layout is not None:
                clear_layout(child_layout)  # Recursively clear nested layouts

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = XChemAlign()
    window.show()

    if CONFIG_FILE.exists():
        window.ui_config.load_yaml(CONFIG_FILE)

    sys.exit(app.exec())
