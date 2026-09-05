#!/usr/bin/env python3
"""Affiche la trajectoire trajkmc.xyz dans une fenêtre OVITO avec
un bouton lecture/pause, un curseur et le numéro de la configuration affichée."""

import argparse
import sys
from pathlib import Path

from PySide6.QtWidgets import QApplication
QApplication(sys.argv)

import ovito
from ovito.io import import_file
from ovito.modifiers import (
    CommonNeighborAnalysisModifier,
    DeleteSelectedModifier,
    ExpressionSelectionModifier,
)
from ovito.vis import Viewport
from ovito.qt_compat import QtCore, QtWidgets

TRAJ_FILE = Path(__file__).with_name("trajkmc.xyz")
PLAYBACK_INTERVAL_MS = 150  # délai entre deux images en lecture automatique

INFO_TEXT = """\
view_trajectory.py — visualiseur de trajectoire KMC (OVITO)

Ce script charge la trajectoire "trajkmc.xyz" (dans le même dossier que
le script) et ouvre une fenêtre OVITO permettant de :

  - visualiser les configurations atomiques successives dans un viewport 3D ;
  - lire la trajectoire automatiquement (bouton "Lecture"/"Pause") ;
  - naviguer manuellement entre les images avec un curseur ;
  - voir le numéro de la configuration affichée (ex. "Configuration : 3 / 42") ;
  - identifier la structure locale de chaque atome via une analyse CNA
    (Common Neighbor Analysis) ;
  - masquer/afficher les atomes de structure FCC (bouton
    "Masquer atomes FCC") pour ne garder que les défauts visibles.

Utilisation :
  python view_trajectory.py            lance la fenêtre de visualisation
  python view_trajectory.py --info     affiche ce message et quitte
"""


class TrajectoryWindow(QtWidgets.QWidget):
    def __init__(self, pipeline):
        super().__init__()
        self.setWindowTitle("Trajectoire KMC")
        self.resize(900, 700)

        self.anim = ovito.scene.anim
        self.first_frame = self.anim.first_frame
        self.last_frame = self.anim.last_frame

        # Analyse des structures locales (CNA) + suppression optionnelle des atomes FCC
        pipeline.modifiers.append(CommonNeighborAnalysisModifier())
        self.select_fcc = ExpressionSelectionModifier(
            expression=f"StructureType=={int(CommonNeighborAnalysisModifier.Type.FCC)}"
        )
        pipeline.modifiers.append(self.select_fcc)
        self.delete_fcc = DeleteSelectedModifier(enabled=False)
        pipeline.modifiers.append(self.delete_fcc)

        vp = Viewport(type=Viewport.Type.Perspective)
        vp.zoom_all()
        self.viewport_widget = vp.create_widget(self)

        self.play_button = QtWidgets.QPushButton("Lecture")
        self.play_button.setCheckable(True)
        self.play_button.toggled.connect(self._toggle_play)

        self.hide_fcc_button = QtWidgets.QPushButton("Masquer atomes FCC")
        self.hide_fcc_button.setCheckable(True)
        self.hide_fcc_button.toggled.connect(self._toggle_hide_fcc)

        self.slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.slider.setMinimum(self.first_frame)
        self.slider.setMaximum(self.last_frame)
        self.slider.setValue(self.anim.current_frame)
        self.slider.valueChanged.connect(self._on_slider_changed)

        self.frame_label = QtWidgets.QLabel()
        self.frame_label.setMinimumWidth(140)
        self._update_frame_label()

        controls = QtWidgets.QHBoxLayout()
        controls.addWidget(self.play_button)
        controls.addWidget(self.hide_fcc_button)
        controls.addWidget(self.slider, stretch=1)
        controls.addWidget(self.frame_label)

        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(self.viewport_widget, stretch=1)
        layout.addLayout(controls)

        self.timer = QtCore.QTimer(self)
        self.timer.setInterval(PLAYBACK_INTERVAL_MS)
        self.timer.timeout.connect(self._advance_frame)

    def _update_frame_label(self):
        self.frame_label.setText(
            f"Configuration : {self.anim.current_frame} / {self.last_frame}"
        )

    def _on_slider_changed(self, value):
        self.anim.current_frame = value
        self._update_frame_label()

    def _advance_frame(self):
        next_frame = self.anim.current_frame + 1
        if next_frame > self.last_frame:
            next_frame = self.first_frame  # reboucle sur la première image
        self.slider.setValue(next_frame)  # déclenche _on_slider_changed

    def _toggle_play(self, checked):
        self.play_button.setText("Pause" if checked else "Lecture")
        if checked:
            self.timer.start()
        else:
            self.timer.stop()

    def _toggle_hide_fcc(self, checked):
        self.delete_fcc.enabled = checked


def parse_args():
    parser = argparse.ArgumentParser(
        description="Visualiseur de trajectoire KMC (OVITO).",
        add_help=True,
    )
    parser.add_argument(
        "--info",
        action="store_true",
        help="afficher une description des fonctionnalités du script et quitter",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    if args.info:
        print(INFO_TEXT)
        return

    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)

    pipeline = import_file(str(TRAJ_FILE))
    pipeline.add_to_scene()

    print(f"Trajectoire chargée : {TRAJ_FILE.name}")
    print(f"Nombre d'images : {pipeline.source.num_frames}")
    print(f"Nombre d'atomes : {pipeline.compute().particles.count}")

    window = TrajectoryWindow(pipeline)
    window.show()

    app.exec()


if __name__ == "__main__":
    main()
