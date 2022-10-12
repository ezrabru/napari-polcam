"""
This module is an example of a barebones QWidget plugin for napari

It implements the Widget specification.
see: https://napari.org/stable/plugins/guides.html?#widgets

Replace code below according to your needs.
"""
import numpy as np
import napari

from vispy.color import Colormap

from typing import TYPE_CHECKING

from magicgui import magic_factory
from qtpy.QtWidgets import QHBoxLayout, QPushButton, QWidget

if TYPE_CHECKING:
    import napari

        
class StokesQWidget(QWidget):
    # your QWidget.__init__ can optionally request the napari viewer instance
    # in one of two ways:
    # 1. use a parameter called `napari_viewer`, as done here
    # 2. use a type annotation of 'napari.viewer.Viewer' for any parameter
    def __init__(self, napari_viewer):
        super().__init__()
        self.viewer = napari_viewer
        
        btn_s0 = QPushButton("Set S0")
        btn_s0.clicked.connect(self._on_click_s0)
        
        btn_s1 = QPushButton("Set S1")
        btn_s1.clicked.connect(self._on_click_s1)
        
        btn_s2 = QPushButton("Set S2")
        btn_s2.clicked.connect(self._on_click_s2)
        
        btn_aolp = QPushButton("Calculate AoLP")
        btn_aolp.clicked.connect(self._on_click_aolp)
        
        btn_dolp = QPushButton("Calculate DoLP")
        btn_dolp.clicked.connect(self._on_click_dolp)
        
        self.setLayout(QHBoxLayout())
        self.layout().addWidget(btn_s0)
        self.layout().addWidget(btn_s1)
        self.layout().addWidget(btn_s2)
        self.layout().addWidget(btn_aolp)
        self.layout().addWidget(btn_dolp)

    def _on_click_s0(self):
        num_layers = 0
        for layer in self.viewer.layers.selection:
            num_layers += 1
        if num_layers == 1:
            # get name of selected layer
            listSelectedLayers = self.viewer.layers.selection
            # save that this layer corresponds to S0
            image = np.random.random((100, 1, 100))
            self.viewer.add_image(image,channel_axis=1,name=["S0"])
            print(f"you have selected S0: {listSelectedLayers}")
        else:
            print(f"You selected {num_layers} layers. Please select 1 layer.")
        
    def _on_click_s1(self):
        num_layers = 0
        for layer in self.viewer.layers.selection:
            num_layers += 1
        if num_layers == 1:
            # get name of selected layer
            listSelectedLayers = self.viewer.layers.selection
            # save that this layer corresponds to S1
            image = np.random.random((100, 1, 100))
            self.viewer.add_image(image,channel_axis=1,name=["S1"])
            print(f"you have selected S1: {listSelectedLayers}")
        else:
            print(f"You selected {num_layers} layers. Please select 1 layer.")
            
    def _on_click_s2(self):
        num_layers = 0
        for layer in self.viewer.layers.selection:
            num_layers += 1
        if num_layers == 1:
            # get name of selected layer
            listSelectedLayers = self.viewer.layers.selection
            # save that this layer corresponds to S2
            image = np.random.random((100, 1, 100))
            self.viewer.add_image(image,channel_axis=1,name=["S2"])
            print(f"you have selected S2: {listSelectedLayers}")
        else:
            print(f"You selected {num_layers} layers. Please select 1 layer.")
            
    def _on_click_aolp(self):
        s1 = self.viewer.layers['S1'].data
        s2 = self.viewer.layers['S2'].data
        AoLP = (1/2)*np.arctan(s2/s1)
        self.viewer.add_image(AoLP,contrast_limits=[-np.pi/2,np.pi/2])
        
    def _on_click_dolp(self):
        s0 = self.viewer.layers['S0'].data
        s1 = self.viewer.layers['S1'].data
        s2 = self.viewer.layers['S2'].data
        DoLP = np.sqrt((s1*s1 + s2*s2)/(s0*s0))
        self.viewer.add_image(DoLP,contrast_limits=[0,1])
        
    def selected_image_layers(self):
        return [
            layer
            for layer in self.viewer.layers.selection
            if isinstance(layer, napari.layers.Image)
        ]


@magic_factory
def example_magic_widget(img_layer: "napari.layers.Image"):
    print(f"you have selected {img_layer}")


# Uses the `autogenerate: true` flag in the plugin manifest
# to indicate it should be wrapped as a magicgui to autogenerate
# a widget.
def example_function_widget(img_layer: "napari.layers.Image"):
    print(f"you have selected {img_layer}")
    

