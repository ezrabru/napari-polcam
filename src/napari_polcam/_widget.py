"""
This module is an example of a barebones QWidget plugin for napari

It implements the Widget specification.
see: https://napari.org/stable/plugins/guides.html?#widgets

Replace code below according to your needs.
"""
import numpy as np
import napari
import pyqtgraph as pg

from matplotlib.colors import hsv_to_rgb

from typing import TYPE_CHECKING

from qtpy.QtWidgets import QSpacerItem, QSizePolicy
from qtpy.QtWidgets import (
    QVBoxLayout,
    QHBoxLayout,
    QWidget,
    QLineEdit,
    QPushButton,
    QLabel,
    QSpinBox,
    )

from qtpy.QtCore import Qt
from superqt import QDoubleRangeSlider

if TYPE_CHECKING:
    import napari


class StokesEstimation(QWidget):
    # your QWidget.__init__ can optionally request the napari viewer instance
    # in one of two ways:
    # 1. use a parameter called `napari_viewer`, as done here
    # 2. use a type annotation of 'napari.viewer.Viewer' for any parameter
    def __init__(self, napari_viewer):
        super().__init__()
        self.viewer = napari_viewer
                
        btn_aolp = QPushButton("Calculate AoLP")
        btn_aolp.clicked.connect(self._on_click_aolp)
        
        btn_dolp = QPushButton("Calculate DoLP")
        btn_dolp.clicked.connect(self._on_click_dolp)

        pixsize_xy = QLineEdit()
        pixsize_xy.setText("1.0")
        pixsize_xy.textChanged.connect(self._on_value_change_pixsize)

        lbl_pixsize_xy = QLabel()
        lbl_pixsize_xy.setText("Voxel size xy:")
        lbl_pixsize_xy.setBuddy(pixsize_xy)
        
        pixsize_z = QLineEdit()
        pixsize_z.setText("1.0")
        pixsize_z.textChanged.connect(self._on_value_change_pixsize)
        
        lbl_pixsize_z = QLabel()
        lbl_pixsize_z.setText("Voxel size z:")
        lbl_pixsize_z.setBuddy(pixsize_z)
        
        self.setLayout(QVBoxLayout())
        self.layout().addWidget(btn_aolp)
        self.layout().addWidget(btn_dolp)
        self.layout().addWidget(lbl_pixsize_xy)
        self.layout().addWidget(pixsize_xy)
        self.layout().addWidget(lbl_pixsize_z)
        self.layout().addWidget(pixsize_z)
        
        self.lineEdit_scale_xy = pixsize_xy
        self.lineEdit_scale_z = pixsize_z

    def _on_value_change_pixsize(self):
        value_xy = float(self.lineEdit_scale_xy.text())
        value_z = float(self.lineEdit_scale_z.text())
        for layer in self.viewer.layers:
            layer.scale = [value_z, value_xy, value_xy]
       
    def _on_click_aolp(self):
        s1 = self.viewer.layers['S1'].data
        s2 = self.viewer.layers['S2'].data
        AoLP = (1/2)*np.arctan(s2/s1)
        self.viewer.add_image(AoLP,contrast_limits=[-np.pi/2,np.pi/2],\
                              scale=(float(self.lineEdit_scale_z.text()),\
                                     float(self.lineEdit_scale_xy.text()),\
                                     float(self.lineEdit_scale_xy.text())))
        
    def _on_click_dolp(self):
        s0 = self.viewer.layers['S0'].data
        s1 = self.viewer.layers['S1'].data
        s2 = self.viewer.layers['S2'].data
        DoLP = np.sqrt((s1*s1 + s2*s2)/(s0*s0))
        self.viewer.add_image(DoLP,contrast_limits=[0,1],\
                              scale=(float(self.lineEdit_scale_z.text()),\
                                     float(self.lineEdit_scale_xy.text()),\
                                     float(self.lineEdit_scale_xy.text())))

class HSVmap(QWidget):
    # your QWidget.__init__ can optionally request the napari viewer instance
    # in one of two ways:
    # 1. use a parameter called `napari_viewer`, as done here
    # 2. use a type annotation of 'napari.viewer.Viewer' for any parameter
    def __init__(self, napari_viewer):
        super().__init__()
        self.viewer = napari_viewer
        
        # check if S0, S1 and S2 are open
        if not (self.check_if_s0_is_loaded() and self.check_if_s1_is_loaded() and self.check_if_s2_is_loaded()):
            lbl_no_s_loaded = QLabel()
            lbl_no_s_loaded.setText("No Stokes parameters are loaded. Load S0, S1 and S2 before using 'Generate HSVmap'.")
        else:
            if not self.check_if_dolp_is_loaded():
                self.calculate_dolp()
            if not self.check_if_aolp_is_loaded():
                self.calculate_aolp()
                
            # A container for S0 histogram plot
            graph_container_s0_hist = QWidget()
            self.s0_hist_widget = pg.GraphicsLayoutWidget()
            self.s0_hist_widget.setBackground(None)
            graph_container_s0_hist.setMaximumHeight(100)
            graph_container_s0_hist.setLayout(QVBoxLayout())
            graph_container_s0_hist.layout().addWidget(self.s0_hist_widget)
            
            # A container for DoLP histogram plot
            graph_container_dolp_hist = QWidget()
            self.dolp_hist_widget = pg.GraphicsLayoutWidget()
            self.dolp_hist_widget.setBackground(None)
            graph_container_dolp_hist.setMaximumHeight(100)
            graph_container_dolp_hist.setLayout(QVBoxLayout())
            graph_container_dolp_hist.layout().addWidget(self.dolp_hist_widget)
            
            # Lateral pixel size
            pixsize_xy = QLineEdit()
            pixsize_xy.setText("1.0")
            pixsize_xy.textChanged.connect(self._on_value_change_pixsize)
            lbl_pixsize_xy = QLabel()
            lbl_pixsize_xy.setText("Voxel size xy:")
            lbl_pixsize_xy.setBuddy(pixsize_xy)
            
            # Axial pixel size
            pixsize_z = QLineEdit()
            pixsize_z.setText("1.0")
            pixsize_z.textChanged.connect(self._on_value_change_pixsize)
            lbl_pixsize_z = QLabel()
            lbl_pixsize_z.setText("Voxel size z:")
            lbl_pixsize_z.setBuddy(pixsize_z)
            
            # DoLP threshold slider
            slider_dolp = QDoubleRangeSlider()
            slider_dolp.setOrientation(Qt.Horizontal)
            slider_dolp.setMinimum(0)
            slider_dolp.setMaximum(1)
            slider_dolp.setValue([0,1])
            slider_dolp.setSingleStep(0.001)
            lbl_slider_dolp = QLabel()
            lbl_slider_dolp.setText("DoLP threshold:")
            lbl_slider_dolp.setAlignment(Qt.AlignCenter)
            lbl_slider_dolp.setBuddy(slider_dolp)
            
            # S0 threshold slider
            slider_s0 = QDoubleRangeSlider()
            slider_s0.setOrientation(Qt.Horizontal)
            slider_s0.setMinimum(0)
            slider_s0.setMaximum(1)
            slider_s0.setValue([0,1])
            slider_s0.setSingleStep(0.01)
            lbl_slider_s0 = QLabel()
            lbl_slider_s0.setText("S0 threshold:")
            lbl_slider_s0.setAlignment(Qt.AlignCenter)
            lbl_slider_s0.setBuddy(slider_s0)
            
            # Button to generate HSVmap on whole dataset
            btn_hsv_map = QPushButton("Calculate HSVmap")
            btn_hsv_map.clicked.connect(self._on_click_hsvmap)
            
            # Assemble all widgets (in order of appearance in gui from top to bottom)
            self.setLayout(QVBoxLayout())
            
            self.layout().addWidget(graph_container_s0_hist)
            self.draw_s0_histogram()
            self.layout().addWidget(lbl_slider_s0)
            self.layout().addWidget(slider_s0)
            
            self.layout().addWidget(graph_container_dolp_hist)
            self.draw_dolp_histogram()
            self.layout().addWidget(lbl_slider_dolp)
            self.layout().addWidget(slider_dolp)
    
            verticalSpacer = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)
            self.layout().addItem(verticalSpacer)
            self.layout().setSpacing(0)
    
            self.layout().addWidget(btn_hsv_map)
    
            verticalSpacer = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)
            self.layout().addItem(verticalSpacer)
            self.layout().setSpacing(0)
            
            self.layout().addWidget(lbl_pixsize_xy)
            self.layout().addWidget(pixsize_xy)
            self.layout().addWidget(lbl_pixsize_z)
            self.layout().addWidget(pixsize_z)
            
            self.lineEdit_scale_xy = pixsize_xy
            self.lineEdit_scale_z = pixsize_z
            
            verticalSpacer = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)
            self.layout().addItem(verticalSpacer)
            self.layout().setSpacing(0)
        
        
    def draw_s0_histogram(self):
        # add a new plot to the s0_hist_widget or empty the old plot
        if not hasattr(self, "plot_s0"):
            self.plot_s0 = self.s0_hist_widget.addPlot()
        else:
            self.plot_s0.clear()
        
        # draw DoLP histogram
        max_s0 = np.max(self.viewer.layers['S0'].data)
        counts, binedges = np.histogram(self.viewer.layers['S0'].data,bins=100,range=(0.0,max_s0))
        self.plot_s0.plot(binedges[0:-1],counts/np.max(counts),name='S0')
        self.plot_s0.setXRange(0, max_s0, padding=0)
        self.plot_s0.setYRange(0, 1, padding=0)
        
    def draw_dolp_histogram(self):
        # add a new plot to the dolp_hist_widget or empty the old plot
        if not hasattr(self, "plot_dolp"):
            self.plot_dolp = self.dolp_hist_widget.addPlot()
        else:
            self.plot_dolp.clear()
            
        # draw DoLP histogram
        counts, binedges = np.histogram(self.viewer.layers['DoLP'].data,bins=100,range=(0.0,1.0))
        self.plot_dolp.plot(binedges[0:-1],counts/np.max(counts),name='DoLP')
        self.plot_dolp.setXRange(0, 1, padding=0)
        self.plot_dolp.setYRange(0, 1, padding=0)
    
    def _on_value_change_pixsize(self):
        value_xy = float(self.lineEdit_scale_xy.text())
        value_z = float(self.lineEdit_scale_z.text())
        for layer in self.viewer.layers:
            layer.scale = [value_z, value_xy, value_xy]
            
    def calculate_aolp(self):
        s1 = self.viewer.layers['S1'].data
        s2 = self.viewer.layers['S2'].data
        AoLP = (1/2)*np.arctan(s2/s1)
        self.viewer.add_image(AoLP,contrast_limits=[-np.pi/2,np.pi/2])
        
    def calculate_dolp(self):
        s0 = self.viewer.layers['S0'].data
        s1 = self.viewer.layers['S1'].data
        s2 = self.viewer.layers['S2'].data
        DoLP = np.sqrt((s1*s1 + s2*s2)/(s0*s0))
        self.viewer.add_image(DoLP,contrast_limits=[0,1])
            
    def check_if_s0_is_loaded(self):
        if 'S0' not in self.viewer.layers: return 0
        else: return 1
        
    def check_if_s1_is_loaded(self):
        if 'S1' not in self.viewer.layers: return 0
        else: return 1
        
    def check_if_s2_is_loaded(self):
        if 'S2' not in self.viewer.layers: return 0
        else: return 1
        
    def check_if_dolp_is_loaded(self):
        if 'DoLP' not in self.viewer.layers: return 0
        else: return 1
        
    def check_if_aolp_is_loaded(self):
        if 'AoLP' not in self.viewer.layers: return 0
        else: return 1

   
    #def _on_click_prepare_hsvmap(self):
    #    # generate 8-bit versions of AoLP and DoLP for faster processing
    
    def _on_click_hsvmap(self):
        
        # get Stokes parameter images (assumed to be already open)
        s0 = self.viewer.layers['S0'].data
        s1 = self.viewer.layers['S1'].data
        s2 = self.viewer.layers['S2'].data
        
        # calculate AoLP (SHOULD CHECK IF ALREADY OPEN)
        AoLP = (1/2)*np.arctan(s2/s1)
        self.viewer.add_image(AoLP,contrast_limits=[-np.pi/2,np.pi/2],\
                              scale=(float(self.lineEdit_scale_z.text()),\
                                     float(self.lineEdit_scale_xy.text()),\
                                     float(self.lineEdit_scale_xy.text())))
        
        # calculate DoLP (SHOULD CHECK IF ALREADY OPEN)
        DoLP = np.sqrt((s1*s1 + s2*s2)/(s0*s0))
        self.viewer.add_image(DoLP,contrast_limits=[0,1],\
                              scale=(float(self.lineEdit_scale_z.text()),\
                                     float(self.lineEdit_scale_xy.text()),\
                                     float(self.lineEdit_scale_xy.text())))
            
        h = self.viewer.layers['AoLP'].data
        s = self.viewer.layers['DoLP'].data
        v = self.viewer.layers['S0'].data
        
        h = (h + np.pi/2)/np.pi; # rescale [-pi/2, pi/2]to [0, 1]
        v = (v - np.min(v))/(np.max(v) - np.min(v)); # rescale intensity to [0 1]
        
        numDim = len(h.shape) # number of dimensions of the dataset
        hsv = np.stack([h,s,v],numDim) # stack
        rgb = hsv_to_rgb(hsv)
        
        hsv_r = rgb[..., 0]*255
        hsv_g = rgb[..., 1]*255
        hsv_b = rgb[..., 2]*255
        
        hsv_r = hsv_r.astype(np.uint8)
        hsv_g = hsv_g.astype(np.uint8)
        hsv_b = hsv_b.astype(np.uint8)
        
        # hide all layers
        #for layer in self.viewer.layers.selection:
        #    layer.selected = False
        
        # add new images for the 3 colour channels
        self.viewer.add_image(hsv_r,contrast_limits=[0,255],colormap="red",blending="additive",\
                              scale=(float(self.lineEdit_scale_z.text()),\
                                     float(self.lineEdit_scale_xy.text()),\
                                     float(self.lineEdit_scale_xy.text())))
        self.viewer.add_image(hsv_g,contrast_limits=[0,255],colormap="green",blending="additive",\
                              scale=(float(self.lineEdit_scale_z.text()),\
                                     float(self.lineEdit_scale_xy.text()),\
                                     float(self.lineEdit_scale_xy.text())))
        self.viewer.add_image(hsv_b,contrast_limits=[0,255],colormap="blue",blending="additive",\
                              scale=(float(self.lineEdit_scale_z.text()),\
                                     float(self.lineEdit_scale_xy.text()),\
                                     float(self.lineEdit_scale_xy.text())))
        self._on_value_change_pixsize # adjust the voxel scaling

class DoLPmap(QWidget):
    # your QWidget.__init__ can optionally request the napari viewer instance
    # in one of two ways:
    # 1. use a parameter called `napari_viewer`, as done here
    # 2. use a type annotation of 'napari.viewer.Viewer' for any parameter
    def __init__(self, napari_viewer):
        super().__init__()
        self.viewer = napari_viewer
        
        pixsize_xy = QLineEdit()
        pixsize_xy.setText("1.0")
        pixsize_xy.textChanged.connect(self._on_value_change_pixsize)

        lbl_pixsize_xy = QLabel()
        lbl_pixsize_xy.setText("Voxel size xy:")
        lbl_pixsize_xy.setBuddy(pixsize_xy)
        
        pixsize_z = QLineEdit()
        pixsize_z.setText("1.0")
        pixsize_z.textChanged.connect(self._on_value_change_pixsize)
        
        lbl_pixsize_z = QLabel()
        lbl_pixsize_z.setText("Voxel size z:")
        lbl_pixsize_z.setBuddy(pixsize_z)
        
        slider_dolp = QDoubleRangeSlider()
        slider_dolp.setOrientation(Qt.Horizontal)
        slider_dolp.setMinimum(0)
        slider_dolp.setMaximum(1)
        slider_dolp.setValue([0,1])
        slider_dolp.setSingleStep(0.001)
        lbl_slider_dolp = QLabel()
        lbl_slider_dolp.setText("DoLP threshold:")
        lbl_slider_dolp.setAlignment(Qt.AlignCenter)
        lbl_slider_dolp.setBuddy(slider_dolp)
        
        slider_s0 = QDoubleRangeSlider()
        slider_s0.setOrientation(Qt.Horizontal)
        slider_s0.setMinimum(0)
        slider_s0.setMaximum(1)
        slider_s0.setValue([0,1])
        slider_s0.setSingleStep(0.01)
        lbl_slider_s0 = QLabel()
        lbl_slider_s0.setText("S0 threshold:")
        lbl_slider_s0.setAlignment(Qt.AlignCenter)
        lbl_slider_s0.setBuddy(slider_s0)
        
        btn_hsv_map = QPushButton("Calculate HSVmap")
        btn_hsv_map.clicked.connect(self._on_click_hsvmap)
        
        self.setLayout(QVBoxLayout())
        
        self.layout().addWidget(btn_hsv_map)
        self.layout().addWidget(lbl_slider_dolp)
        self.layout().addWidget(slider_dolp)
        self.layout().addWidget(lbl_slider_s0)
        self.layout().addWidget(slider_s0)
        
        self.layout().addWidget(lbl_pixsize_xy)
        self.layout().addWidget(pixsize_xy)
        self.layout().addWidget(lbl_pixsize_z)
        self.layout().addWidget(pixsize_z)
        
        self.lineEdit_scale_xy = pixsize_xy
        self.lineEdit_scale_z = pixsize_z

    def _on_value_change_pixsize(self):
        value_xy = float(self.lineEdit_scale_xy.text())
        value_z = float(self.lineEdit_scale_z.text())
        for layer in self.viewer.layers:
            layer.scale = [value_z, value_xy, value_xy]
       
    #def _on_click_prepare_hsvmap(self):
    #    # generate 8-bit versions of AoLP and DoLP for faster processing
    
    def _on_click_hsvmap(self):
        
        # get Stokes parameter images (assumed to be already open)
        s0 = self.viewer.layers['S0'].data
        s1 = self.viewer.layers['S1'].data
        s2 = self.viewer.layers['S2'].data
        
        # calculate AoLP (SHOULD CHECK IF ALREADY OPEN)
        AoLP = (1/2)*np.arctan(s2/s1)
        self.viewer.add_image(AoLP,contrast_limits=[-np.pi/2,np.pi/2],\
                              scale=(float(self.lineEdit_scale_z.text()),\
                                     float(self.lineEdit_scale_xy.text()),\
                                     float(self.lineEdit_scale_xy.text())))
        
        # calculate DoLP (SHOULD CHECK IF ALREADY OPEN)
        DoLP = np.sqrt((s1*s1 + s2*s2)/(s0*s0))
        self.viewer.add_image(DoLP,contrast_limits=[0,1],\
                              scale=(float(self.lineEdit_scale_z.text()),\
                                     float(self.lineEdit_scale_xy.text()),\
                                     float(self.lineEdit_scale_xy.text())))
            
        h = self.viewer.layers['AoLP'].data
        s = self.viewer.layers['DoLP'].data
        v = self.viewer.layers['S0'].data
        
        h = (h + np.pi/2)/np.pi; # rescale [-pi/2, pi/2]to [0, 1]
        v = (v - np.min(v))/(np.max(v) - np.min(v)); # rescale intensity to [0 1]
        
        numDim = len(h.shape) # number of dimensions of the dataset
        hsv = np.stack([h,s,v],numDim) # stack
        rgb = hsv_to_rgb(hsv)
        
        hsv_r = rgb[..., 0]*255
        hsv_g = rgb[..., 1]*255
        hsv_b = rgb[..., 2]*255
        
        hsv_r = hsv_r.astype(np.uint8)
        hsv_g = hsv_g.astype(np.uint8)
        hsv_b = hsv_b.astype(np.uint8)
        
        # hide all layers
        #for layer in self.viewer.layers.selection:
        #    layer.selected = False
        
        # add new images for the 3 colour channels
        self.viewer.add_image(hsv_r,contrast_limits=[0,255],colormap="red",blending="additive",\
                              scale=(float(self.lineEdit_scale_z.text()),\
                                     float(self.lineEdit_scale_xy.text()),\
                                     float(self.lineEdit_scale_xy.text())))
        self.viewer.add_image(hsv_g,contrast_limits=[0,255],colormap="green",blending="additive",\
                              scale=(float(self.lineEdit_scale_z.text()),\
                                     float(self.lineEdit_scale_xy.text()),\
                                     float(self.lineEdit_scale_xy.text())))
        self.viewer.add_image(hsv_b,contrast_limits=[0,255],colormap="blue",blending="additive",\
                              scale=(float(self.lineEdit_scale_z.text()),\
                                     float(self.lineEdit_scale_xy.text()),\
                                     float(self.lineEdit_scale_xy.text())))
        self._on_value_change_pixsize # adjust the voxel scaling
        
    def selected_image_layers(self):
        return [
            layer
            for layer in self.viewer.layers.selection
            if isinstance(layer, napari.layers.Image)
        ]

