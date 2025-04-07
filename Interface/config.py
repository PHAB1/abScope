import sys
from config import *
from PyQt6.QtWidgets import (
    QMainWindow, QApplication, QFileDialog, QToolBar, QStatusBar, QMessageBox, QPushButton, QWidget, QVBoxLayout, QLabel, QCheckBox, QDialog
)
from PyQt6.QtGui import QAction, QIcon, QPalette, QColor
from PyQt6.QtCore import Qt

class ConfigWindow(QDialog):  # QDialog é ideal para pop-ups
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Configurações")
        self.setGeometry(200, 200, 300, 200)

        layout = QVBoxLayout()

        self.label = QLabel("Opções:")
        layout.addWidget(self.label)

        self.option1 = QCheckBox("Ativar Modo Escuro")
        layout.addWidget(self.option1)

        self.option2 = QCheckBox("Habilitar Notificações")
        layout.addWidget(self.option2)

        self.close_button = QPushButton("Fechar")
        self.close_button.clicked.connect(self.close)
        layout.addWidget(self.close_button)

        self.setLayout(layout)
