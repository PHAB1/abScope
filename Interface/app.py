import sys
import os
from config import *
from PyQt6.QtWidgets import (
    QMainWindow, QApplication, QFileDialog, QToolBar, QStatusBar, QLabel, 
    QWidget, QVBoxLayout, QListWidget, QHBoxLayout, QListWidgetItem,
    QPushButton, QSizePolicy
)
from PyQt6.QtGui import QAction, QIcon, QPalette, QColor, QLinearGradient, QPainter
from PyQt6.QtCore import Qt, QSize, QPoint

class CustomListWidget(QListWidget):
    def __init__(self):
        super().__init__()
        self.file_paths = []
        self.setMaximumWidth(300)
        self.setStyleSheet("""
            QListWidget {
                background-color: #282a36;
                border: none;
                border-radius: 8px;
                padding: 5px;
                margin: 10px;
                color: #f8f8f2;
                font-size: 14px;
                outline: none;
            }
            QListWidget::item {
                border-radius: 5px;
                margin: 2px;
            }
            QListWidget::item:hover {
                background-color: #343744;
            }
            QListWidget::item:selected {
                background-color: #44475a;
                border: 1px solid #6272a4;
            }
        """)

    def add_files(self, files):
        for file_path in files:
            if file_path not in self.file_paths:
                self._add_file_item(file_path)
                self.file_paths.append(file_path)

    def _add_file_item(self, file_path):
        item = QListWidgetItem()
        item.setSizeHint(QSize(200, 50))
        
        widget = QWidget()
        widget.setStyleSheet("background: transparent;")
        layout = QHBoxLayout()
        layout.setContentsMargins(10, 5, 10, 5)
        layout.setSpacing(10)
        
        label = QLabel(os.path.basename(file_path))
        label.setToolTip(file_path)
        label.setStyleSheet("""
            QLabel {
                color: #f8f8f2;
                font-size: 13px;
                background: transparent;
            }
        """)
        
        btn = QPushButton("×")
        btn.setFixedSize(24, 24)
        btn.setStyleSheet("""
            QPushButton {
                background-color: #ff5555;
                color: white;
                border-radius: 4px;
                font-weight: bold;
                border: none;
                font-size: 16px;
                padding-bottom: 2px;
            }
            QPushButton:hover {
                background-color: #ff6e6e;
                border-radius: 4px;
            }
            QPushButton:pressed {
                background-color: #cc4444;
            }
        """)
        btn.clicked.connect(lambda _, i=item, fp=file_path: self.remove_item(i, fp))
        
        layout.addWidget(label, 1)
        layout.addWidget(btn)
        widget.setLayout(layout)
        
        self.addItem(item)
        self.setItemWidget(item, widget)

    def remove_item(self, item, file_path):
        row = self.row(item)
        self.file_paths.remove(file_path)
        self.takeItem(row)

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Modern File Manager")
        self.setGeometry(100, 100, 1000, 700)
        self.apply_modern_theme()
        self.setAcceptDrops(True)

        # Configuração do layout principal
        main_layout = QHBoxLayout()
        main_layout.setSpacing(15)
        main_layout.setContentsMargins(15, 15, 15, 15)
        
        self.file_list_widget = CustomListWidget()
        
        main_content = QWidget()
        main_content.setStyleSheet("""
            background-color: #44475a;
            border-radius: 10px;
        """)
        
        main_layout.addWidget(self.file_list_widget)
        main_layout.addWidget(main_content, 1)
        
        container = QWidget()
        container.setLayout(main_layout)
        self.setCentralWidget(container)
        
        # Configuração da toolbar
        self.toolbar = QToolBar("Main Toolbar")
        self.toolbar.setMovable(False)
        self.toolbar.setStyleSheet("""
            QToolBar {
                background: qlineargradient(x1:0, y1:0, x1:1, y1:0,
                    stop:0 #282a36, stop:1 #44475a);
                padding: 8px;
                spacing: 10px;
                border-bottom: 1px solid #6272a4;
            }
            QToolButton {
                background: #6272a4;
                border-radius: 5px;
                padding: 8px;
                color: #f8f8f2;
            }
            QToolButton:hover {
                background: #7282b4;
            }
            QToolButton:pressed {
                background: #525284;
            }
        """)
        self.addToolBar(self.toolbar)
        
        # Ações da toolbar
        self.create_action("Abrir", "document-open", self.open_files_dialog, "Ctrl+O")
        self.create_action("Configurações", "preferences-system", self.config_file, "Ctrl+P")
        
        # Status bar
        self.statusBar().setStyleSheet("""
            QStatusBar {
                background-color: #282a36;
                color: #6272a4;
                border-top: 1px solid #6272a4;
                padding-left: 8px;
            }
        """)
        self.statusBar().showMessage("Pronto")

    def create_action(self, text, icon_name, callback, shortcut=None):
        action = QAction(QIcon.fromTheme(icon_name), text, self)
        action.triggered.connect(callback)
        if shortcut:
            action.setShortcut(shortcut)
        self.toolbar.addAction(action)
        return action

    def apply_modern_theme(self):
        palette = QPalette()
        palette.setColor(QPalette.ColorRole.Window, QColor("#282a36"))
        palette.setColor(QPalette.ColorRole.WindowText, QColor("#f8f8f2"))
        palette.setColor(QPalette.ColorRole.Base, QColor("#44475a"))
        palette.setColor(QPalette.ColorRole.Text, QColor("#f8f8f2"))
        palette.setColor(QPalette.ColorRole.Button, QColor("#6272a4"))
        palette.setColor(QPalette.ColorRole.ButtonText, QColor("#f8f8f2"))
        palette.setColor(QPalette.ColorRole.Highlight, QColor("#bd93f9"))
        palette.setColor(QPalette.ColorRole.HighlightedText, QColor("#282a36"))
        self.setPalette(palette)

        self.setStyleSheet("""
            QMainWindow {
                background-color: #282a36;
            }
            QScrollBar:vertical {
                width: 12px;
                background: #343746;
                border-radius: 6px;
            }
            QScrollBar::handle:vertical {
                background: #6272a4;
                min-height: 20px;
                border-radius: 6px;
            }
            QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {
                height: 0px;
            }
        """)

    def open_files_dialog(self):
        files, _ = QFileDialog.getOpenFileNames(self, "Selecionar Arquivos")
        if files:
            self.file_list_widget.add_files(files)
            self.statusBar().showMessage(f"{len(files)} arquivo(s) adicionado(s)")

    def config_file(self):
        self.config_window = ConfigWindow()
        self.config_window.exec()

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls():
            event.acceptProposedAction()

    def dropEvent(self, event):
        files = [url.toLocalFile() for url in event.mimeData().urls()]
        self.file_list_widget.add_files(files)
        self.statusBar().showMessage(f"{len(files)} arquivo(s) adicionado(s) via drag & drop")

app = QApplication(sys.argv)
app.setStyle("Fusion")
window = MainWindow()
window.show()
app.exec()