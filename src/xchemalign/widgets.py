
from PyQt5.QtGui import QFont
from PyQt5.QtWidgets import QLineEdit, QTextEdit, QComboBox

MONO_FONT = QFont("Courier New")
MONO_FONT.setStyleHint(QFont.Monospace)

class ModStyleMixin:
    def style_modification(self):
        self.setStyleSheet("font-weight:bold; color:orange")

    def style_reset(self):
        self.setStyleSheet("")

class LineEdit(QLineEdit,ModStyleMixin):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setFont(MONO_FONT)
        self.style_reset()
        self.textEdited.connect(self.style_modification)

class TextEdit(QTextEdit,ModStyleMixin):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setFont(MONO_FONT)
        self.style_reset()        
        self.setMinimumHeight(100)
        self.textChanged.connect(self.style_modification)

    def setPlainText(self, *args, **kwargs):
        super().setPlainText(*args, **kwargs)
        self.style_reset()

class ComboBox(QComboBox,ModStyleMixin):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setFont(MONO_FONT)
        self.style_reset()
        self.starting_index = int(self.currentIndex())

    def connect(self):
        self.currentIndexChanged.connect(self.style_modification)

    def style_modification(self):
        if self.currentIndex() != self.starting_index:
            self.setStyleSheet("font-weight:bold; color:orange")
        else:
            self.style_reset()

