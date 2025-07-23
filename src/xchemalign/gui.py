import sys

import subprocess
import signal

signal.signal(signal.SIGINT, signal.SIG_DFL)

from datetime import datetime

from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont

from PyQt5.QtWidgets import (
    QGroupBox,
    QApplication,
    # QCheckBox,
    QComboBox,
    # QDateEdit,
    # QDateTimeEdit,
    # QDial,
    # QDoubleSpinBox,
    # QFontComboBox,
    QLabel,
    # QLCDNumber,
    QLineEdit,
    # QMainWindow,
    # QProgressBar,
    QPushButton,
    # QRadioButton,
    # QSlider,
    # QSpinBox,
    # QTimeEdit,
    QHBoxLayout,
    QVBoxLayout,
    # QFrame,
    QTextEdit,
    QWidget,
    QScrollArea,
)

import pandas as pd
import mrich
from pathlib import Path
from shutil import copy2

from .widgets import *

DIRECTORY = Path(".").resolve()
CONFIG_FILE = Path("config.yaml").resolve()
ASSEMBLIES_FILE = Path("assemblies.yaml").resolve()

CONFIG_BACKUP = None
ASSEMBLIES_BACKUP = None

class XChemAlign(QWidget):

    def __init__(self):
        super().__init__()
        self.setup_window()
        self.setup_layout()
        self.ui_config = ConfigUI(self, self.layout_config)
        self.ui_assemblies = AssembliesUI(self, self.layout_assemblies)
        self.ui_crystalforms = CrystalformsUI(self, self.layout_crystalforms)

    def setup_window(self, width=1200, height=800):
        self.setWindowTitle("XChemAlign GUI")
        self.resize(width, height)

    def setup_layout(self):

        self.layout = QVBoxLayout()

        ### HEADER
        self.layout_header = QVBoxLayout()
        label = QLabel(f"{DIRECTORY}")
        label.setFont(MONO_FONT)
        label.setStyleSheet("background-color:black; color:white; padding:8px")
        self.layout_header.addWidget(label)

        button_layout = QHBoxLayout()
        save = QPushButton("Save Changes")
        button_layout.addWidget(save)
        revert = QPushButton("Revert to backup")
        revert.clicked.connect(revert_state)
        button_layout.addWidget(revert)
        collate = QPushButton("Run Collator")
        collate.clicked.connect(run_collator)
        button_layout.addWidget(collate)
        align = QPushButton("Run Aligner")
        align.clicked.connect(run_aligner)
        button_layout.addWidget(align)
        self.layout_header.addLayout(button_layout)

        self.layout.addLayout(self.layout_header)
        
        ### MAIN LAYOUT
        self.layout_bottom = QHBoxLayout()
        self.layout.addLayout(self.layout_bottom)

        # sublayouts
        self.layout_left = QVBoxLayout()
        self.layout_middle = QVBoxLayout()
        # self.layout_right = QVBoxLayout()
        self.layout_bottom.addLayout(self.layout_left)
        self.layout_bottom.addLayout(self.layout_middle)
        # self.layout_bottom.addLayout(self.layout_right)

        # config layout
        self.layout_config = QVBoxLayout()
        self.layout_left.addLayout(self.layout_config)
        label = QLabel("CONFIG")
        font = label.font()
        font.setPointSize(14)
        font.setBold(True)
        label.setAlignment(Qt.AlignCenter)
        label.setFont(font)
        self.layout_config.addWidget(label)
                
        # assemblies layout
        self.layout_assemblies = QVBoxLayout()
        self.layout_middle.addLayout(self.layout_assemblies)
        label = QLabel("ASSEMBLIES")
        font = label.font()
        font.setPointSize(14)
        font.setBold(True)
        label.setAlignment(Qt.AlignCenter)
        label.setFont(font)
        self.layout_assemblies.addWidget(label)

        # crystalforms layout
        self.layout_crystalforms = QVBoxLayout()
        self.layout_middle.addLayout(self.layout_crystalforms)
        label = QLabel("CRYSTALFORMS")
        font = label.font()
        font.setPointSize(14)
        font.setBold(True)
        label.setAlignment(Qt.AlignCenter)
        label.setFont(font)
        self.layout_crystalforms.addWidget(label)

        self.setLayout(self.layout)


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
        self.widget_target_name = LineEdit()
        l = QHBoxLayout()
        l.addWidget(self.label_target_name)
        l.addWidget(self.widget_target_name)
        self.layout.addLayout(l)

        # target name
        self.label_base_dir = QLabel("Base Directory:")
        self.widget_base_dir = LineEdit("/")
        l = QHBoxLayout()
        l.addWidget(self.label_base_dir)
        l.addWidget(self.widget_base_dir)
        self.widget_base_dir.setFont(MONO_FONT)
        self.layout.addLayout(l)

        # inputs
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
                soakdb=d.get("soakdb", ""),
                panddas_event_files=d.get("panddas_event_files", ""),
                exclude=d.get("exclude", ""),
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

        self.button_add = QPushButton("Add input")
        self.button_add.clicked.connect(self.add_input)
        
        l = QHBoxLayout()
        l.addWidget(self.label_title)
        l.addWidget(self.button_add)
        self.layout.addLayout(l)

        self.layout_inputs = QVBoxLayout()
        self.layout.addLayout(self.layout_inputs)

        # self.layout_buttons = QHBoxLayout()

        # b1 = QPushButton("SAVE config.yaml")
        # b1.setStyleSheet("color:green")

        # def b1_action():
        #     raise NotImplementedError

        # b1.clicked.connect(b1_action)
        # self.layout_buttons.addWidget(b1)

        # b2 = QPushButton("REVERT config.yaml")
        # b2.setStyleSheet("color:red")

        # def b2_action():
        #     self.parent.load_yaml(CONFIG_FILE)

        # b2.clicked.connect(b2_action)
        # self.layout_buttons.addWidget(b2)

        # self.layout.addLayout(self.layout_buttons)

        # ###

        self.update_widgets()

    def add_input(
        self,
        *args,
        dir="",
        type="model_building",
        code_prefix="",
        code_prefix_tooltip="",
        soakdb="",
        panddas_event_files="",
        exclude="",
    ):

        df = pd.DataFrame(
            [
                dict(
                    dir=dir,
                    type=type,
                    code_prefix=code_prefix,
                    code_prefix_tooltip=code_prefix_tooltip,
                    soakdb=soakdb,
                    panddas_event_files=panddas_event_files,
                    exclude=exclude,
                )
            ]
        )

        self.inputs = pd.concat([self.inputs, df], ignore_index=True)

        self.update_widgets()

    def remove_input(self, *args, index=None):
        self.inputs = self.inputs.drop(index)
        self.update_widgets()

    def update_widgets(self):

        clear_layout(self.layout_inputs)

        container_widget = QWidget()
        inner_layout = QVBoxLayout(container_widget)

        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setWidget(container_widget)

        self.layout_inputs.addWidget(scroll_area)

        for i, row in self.inputs.iterrows():

            groupbox = QGroupBox()
            layout = QVBoxLayout()
            groupbox.setLayout(layout)
            inner_layout.addWidget(groupbox)

            # title
            l = QHBoxLayout()
            label = QLabel(f"INPUT {i+1}")
            font = label.font()
            # font.setPointSize(14)
            font.setBold(True)
            # label.setAlignment(Qt.AlignCenter)
            label.setFont(font)
            l.addWidget(label)

            # remove
            button = QPushButton("remove")
            button.clicked.connect(lambda x, idx=i: self.remove_input(index=idx))
            button.setStyleSheet("color: red;")
            l.addWidget(button)

            layout.addLayout(l)

            # dir
            label = QLabel("dir")
            widget = LineEdit(row["dir"])
            l = QHBoxLayout()
            l.addWidget(label)
            l.addWidget(widget)
            layout.addLayout(l)

            # type
            label = QLabel("type")
            widget = ComboBox()
            widget.addItems(["model_building", "manual"])
            index = widget.findText(row["type"])
            if index != -1:
                widget.setCurrentIndex(index)
                widget.starting_index = int(index)
            widget.connect()

            l = QHBoxLayout()
            l.addWidget(label)
            l.addWidget(widget)
            layout.addLayout(l)

            # code_prefix
            label = QLabel("code_prefix")
            widget = LineEdit(row["code_prefix"])
            l = QHBoxLayout()
            l.addWidget(label)
            l.addWidget(widget)
            layout.addLayout(l)

            # code_prefix
            label = QLabel("code_prefix_tooltip")
            widget = LineEdit(row["code_prefix_tooltip"])
            l = QHBoxLayout()
            l.addWidget(label)
            l.addWidget(widget)
            layout.addLayout(l)

            if row["type"] == "model_building":
                # panddas_event_files
                label = QLabel("panddas_event_files")
                widget = TextEdit()
                widget.setPlainText("\n".join(row["panddas_event_files"]))
                layout.addWidget(label)
                layout.addWidget(widget)

            # panddas_event_files
            label = QLabel("exclude")
            widget = TextEdit()
            widget.setPlainText("\n".join(row["exclude"]))
            layout.addWidget(label)
            layout.addWidget(widget)

class AssembliesUI(ChildUI):

    def clear(self):
        self.assemblies = pd.DataFrame(columns=["assembly", "reference", "biomol", "chains"])

    def setup(self):
        self.clear()

        self.label_title = QLabel("Assemblies")
        font = self.label_title.font()
        font.setBold(True)
        self.label_title.setFont(font)

        self.button_add = QPushButton("Add assembly")
        self.button_add.clicked.connect(self.add_assembly)
        
        l = QHBoxLayout()
        l.addWidget(self.label_title)
        l.addWidget(self.button_add)
        self.layout.addLayout(l)

        self.layout_assemblies = QVBoxLayout()
        self.layout.addLayout(self.layout_assemblies)

        self.update_widgets()

    def add_assembly(self, *args, assembly="", reference="", biomol="", chains=""):
        df = pd.DataFrame([dict(assembly=assembly, reference=reference, biomol=biomol, chains=chains)])
        self.assemblies = pd.concat([self.assemblies, df], ignore_index=True)
        self.update_widgets()

    def remove_assembly(self, *args, index=None):
        self.assemblies = self.assemblies.drop(index)
        self.update_widgets()


    def update_widgets(self):
        
        clear_layout(self.layout_assemblies)

        container_widget = QWidget()
        inner_layout = QVBoxLayout(container_widget)

        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setWidget(container_widget)

        self.layout_assemblies.addWidget(scroll_area)

        for i, row in self.assemblies.iterrows():

            groupbox = QGroupBox()
            layout = QVBoxLayout()
            groupbox.setLayout(layout)
            inner_layout.addWidget(groupbox)

            # title
            l = QHBoxLayout()
            label = QLabel(f"ASSEMBLY {i+1}")
            font = label.font()
            font.setBold(True)
            label.setFont(font)
            l.addWidget(label)

            # remove
            button = QPushButton("remove")
            button.clicked.connect(lambda x, idx=i: self.remove_assembly(index=idx))
            button.setStyleSheet("color: red;")
            l.addWidget(button)

            layout.addLayout(l)

            # assembly name
            l = QHBoxLayout()
            label = QLabel("assembly")
            widget = LineEdit(row["assembly"])
            l.addWidget(label)
            l.addWidget(widget)
            layout.addLayout(l)

            # reference
            l = QHBoxLayout()
            label = QLabel("reference")
            widget = LineEdit(row["reference"])
            l.addWidget(label)
            l.addWidget(widget)
            layout.addLayout(l)

            # biomol
            l = QHBoxLayout()
            label = QLabel("biomol")
            widget = LineEdit(row["biomol"])
            l.addWidget(label)
            l.addWidget(widget)
            layout.addLayout(l)

            # chains
            l = QHBoxLayout()
            label = QLabel("chains")
            widget = LineEdit(row["chains"])
            l.addWidget(label)
            l.addWidget(widget)
            layout.addLayout(l)

    def load_yaml(self, file):

        import yaml

        with open(file, "r") as stream:
            config = yaml.safe_load(stream)

        self.clear()

        for k,d in config["assemblies"].items():
            self.add_assembly(
                assembly=k,
                reference=d["reference"],
                biomol=d["biomol"],
                chains=d["chains"],
            )

class CrystalformsUI(ChildUI):

    def clear(self):
        self.crystalforms = pd.DataFrame(columns=["crystalform", "reference", "assemblies"])

    def setup(self):
        self.clear()

        self.label_title = QLabel("Crystalforms")
        font = self.label_title.font()
        font.setBold(True)
        self.label_title.setFont(font)

        self.button_add = QPushButton("Add crystalform")
        self.button_add.clicked.connect(self.add_crystalform)
        
        l = QHBoxLayout()
        l.addWidget(self.label_title)
        l.addWidget(self.button_add)
        self.layout.addLayout(l)

        self.layout_crystalforms = QVBoxLayout()
        self.layout.addLayout(self.layout_crystalforms)

        # self.update_widgets()

    def add_crystalform(self, *args, crystalform="", reference="", assemblies=""):
        df = pd.DataFrame([dict(crystalform=crystalform, reference=reference, assemblies=assemblies)])
        self.crystalforms = pd.concat([self.crystalforms, df], ignore_index=True)
        self.update_widgets()

    def remove_crystalform(self, *args, index=None):
        self.crystalforms = self.crystalforms.drop(index)
        self.update_widgets()


    def update_widgets(self):
        
        clear_layout(self.layout_crystalforms)

        container_widget = QWidget()
        inner_layout = QVBoxLayout(container_widget)

        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setWidget(container_widget)

        self.layout_crystalforms.addWidget(scroll_area)

        for i, row in self.crystalforms.iterrows():

            groupbox = QGroupBox()
            layout = QVBoxLayout()
            groupbox.setLayout(layout)
            inner_layout.addWidget(groupbox)

            # title
            l = QHBoxLayout()
            label = QLabel(f"CRYSTALFORM {i+1}")
            font = label.font()
            font.setBold(True)
            label.setFont(font)
            l.addWidget(label)

            # remove
            button = QPushButton("remove")
            button.clicked.connect(lambda x, idx=i: self.remove_assembly(index=idx))
            button.setStyleSheet("color: red;")
            l.addWidget(button)

            layout.addLayout(l)

            # crystalform name
            l = QHBoxLayout()
            label = QLabel("crystalform")
            widget = QLineEdit(row["crystalform"])
            l.addWidget(label)
            l.addWidget(widget)
            layout.addLayout(l)

            # reference
            l = QHBoxLayout()
            label = QLabel("reference")
            widget = QLineEdit(row["reference"])
            l.addWidget(label)
            l.addWidget(widget)
            layout.addLayout(l)

            # panddas_event_files
            label = QLabel("assemblies")
            widget = QTextEdit()
            widget.setFont(MONO_FONT)
            # widget.setPlainText("\n".join(row["assemblies"]))
            widget.setPlainText(row["assemblies"])
            layout.addWidget(label)
            layout.addWidget(widget)

    def load_yaml(self, file):

        import yaml

        with open(file, "r") as stream:
            config = yaml.safe_load(stream)

        self.clear()

        for k,d in config["crystalforms"].items():
            self.add_crystalform(
                crystalform=k,
                reference=d["reference"],
                assemblies=yaml.dump(d["assemblies"]),
            )

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

def backup_yaml(file):
    timestamp_str = datetime.now().strftime("%Y-%m-%d.%H-%M-%S")
    new_file = DIRECTORY / file.name.replace(".yaml", f".{timestamp_str}.yaml")
    mrich.writing(new_file)
    copy2(file, new_file)
    return new_file

def revert_state():
    return load_yamls(from_backup=True)

def load_yamls(from_backup=False):

    global CONFIG_BACKUP
    global ASSEMBLIES_BACKUP

    if from_backup: 
        window.ui_config.load_yaml(CONFIG_BACKUP)
        window.ui_assemblies.load_yaml(ASSEMBLIES_BACKUP)
        window.ui_crystalforms.load_yaml(ASSEMBLIES_BACKUP)
        return

    if CONFIG_FILE.exists():
        CONFIG_BACKUP = backup_yaml(CONFIG_FILE)
        window.ui_config.load_yaml(CONFIG_FILE)
    else:
        window.ui_config.ui_inputs.add_input()

    if ASSEMBLIES_FILE.exists():
        ASSEMBLIES_BACKUP = backup_yaml(ASSEMBLIES_FILE)
        window.ui_assemblies.load_yaml(ASSEMBLIES_FILE)
        window.ui_crystalforms.load_yaml(ASSEMBLIES_FILE)
    else:
        window.ui_assemblies.add_assembly()
        window.ui_crystalforms.add_crystalform()

def run_collator():
    subprocess.run(["python", "-m", "xchemalign.collator"])

def run_aligner():
    subprocess.run(["python", "-m", "xchemalign.aligner"])

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = XChemAlign()
    window.show()

    load_yamls()

    sys.exit(app.exec())
