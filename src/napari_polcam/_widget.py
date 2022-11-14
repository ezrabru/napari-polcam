"""
This module is an example of a barebones QWidget plugin for napari

It implements the Widget specification.
see: https://napari.org/stable/plugins/guides.html?#widgets

"""

import numpy as np
from napari.utils.notifications import show_info
from napari.experimental import link_layers
#from napari.utils import progress
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
    QComboBox,  
    QGroupBox,
    QCheckBox
    )

from qtpy.QtCore import Qt
#from superqt import QDoubleRangeSlider

if TYPE_CHECKING:
    import napari

from ._functions import PolarisationCameraImage


class StokesEstimation(QWidget):
    
    def __init__(self, napari_viewer):
        super().__init__()
        self.viewer = napari_viewer
        self.offset = 0.0 # initialise background to be zero
        
        # =====================================================================
        # Settings group 
        # =====================================================================
        
        # layout of the "polariser unit" gui input ============================
        polariser_unit = QWidget()
        polariser_unit.setLayout(QHBoxLayout())
        polariser_unit.setMaximumHeight(40)
        lbl_unit = QLabel("Polariser unit: ")
        polariser_unit.layout().addWidget(lbl_unit)
        dropdown_unit = QComboBox()
        dropdown_unit.addItem("[-45 0; 90 45]")
        dropdown_unit.addItem("[0 -45; 45 90]")
        dropdown_unit.addItem("[90 45; -45 0]")
        dropdown_unit.addItem("[45 90; 0 -45]")
        polariser_unit.layout().addWidget(dropdown_unit)
        polariser_unit.layout().setSpacing(0)
        self.dropdown_unit = dropdown_unit

        # layout of the "background" gui input ================================
        bkgnd_container = QWidget()
        bkgnd_container.setLayout(QHBoxLayout())
        bkgnd_container.setMaximumHeight(40)
        lbl_bkgnd = QLabel("Background:")
        bkgnd_container.layout().addWidget(lbl_bkgnd)
        lineedit_bkgnd = QLineEdit()
        lineedit_bkgnd.setText(str(self.offset))
        lineedit_bkgnd.textChanged.connect(self._on_value_change_bkgnd)
        bkgnd_container.layout().addWidget(lineedit_bkgnd)
        self.lineedit_bkgnd = lineedit_bkgnd

        # layout of the "interpolation method" gui input ======================
        method_choice = QWidget()
        method_choice.setLayout(QHBoxLayout())
        lbl_method = QLabel("Method: ")
        method_choice.layout().addWidget(lbl_method)
        dropdown_method = QComboBox()
        dropdown_method.addItem("Cubic spline interpolation")
        #dropdown_method.addItem("Fourier")
        dropdown_method.addItem("None")
        method_choice.layout().addWidget(dropdown_method)
        self.dropdown_method = dropdown_method
        
        # checkbox for showing intermediate results ===========================
        checkbox_show_intermediate_results = QCheckBox("Show intermediate results as new layers")
        checkbox_show_intermediate_results.setCheckState(False)
        self.checkbox_show_intermediate_results = checkbox_show_intermediate_results
        
        # group all settings gui elements in a box ============================
        settingsGroupBox = QGroupBox("Settings")
        settings_box = QVBoxLayout()
        settings_box.addWidget(polariser_unit)        
        settings_box.addWidget(method_choice)
        settings_box.addWidget(bkgnd_container)
        settings_box.addWidget(checkbox_show_intermediate_results)
        settings_box.addStretch(1)
        settings_box.setSpacing(0)
        settingsGroupBox.setLayout(settings_box)
        
        
        # =====================================================================
        # Basic processing group 
        # =====================================================================
        
        btn_channels = QPushButton("Calculate channels")
        btn_channels.clicked.connect(self._on_click_channels)
        
        btn_stokes = QPushButton("Calculate Stokes parameters")
        btn_stokes.clicked.connect(self._on_click_stokes_add_to_viewer)
        
        btn_quadview = QPushButton("Calculate Quadview")
        btn_quadview.clicked.connect(self._on_click_quadview)
        
        # A container for AoLP and DoLP buttons ===============================
        btn_aolp_dolp = QWidget()
        btn_aolp_dolp.setLayout(QHBoxLayout())
        btn_aolp = QPushButton("Calculate AoLP")
        btn_dolp = QPushButton("Calculate DoLP")
        btn_aolp.clicked.connect(self._on_click_aolp)        
        btn_dolp.clicked.connect(self._on_click_dolp)
        btn_aolp_dolp.layout().addWidget(btn_aolp)
        btn_aolp_dolp.layout().addWidget(btn_dolp)
        
        # group all basic processing gui elements in a box ====================
        basicProcessingGroupBox = QGroupBox("Basic processing")
        basic_processing_box = QVBoxLayout()
        basic_processing_box.addWidget(btn_quadview)
        basic_processing_box.addWidget(btn_channels)        
        basic_processing_box.addWidget(btn_stokes)
        basic_processing_box.addWidget(btn_aolp_dolp)
        basic_processing_box.addStretch(1)
        basicProcessingGroupBox.setLayout(basic_processing_box)
        
        
        # =====================================================================
        # AoLP/DoLP/S0 HSV colourmap rendering group 
        # =====================================================================
        
        colmap_choice_container = QWidget()
        colmap_choice_container.setLayout(QHBoxLayout())
        lbl_colmap_choice = QLabel("Colourmap: ")
        colmap_choice_container.layout().addWidget(lbl_colmap_choice)
        
        dropdown_colmap = QComboBox()
        dropdown_colmap.addItem("HSVmap")
        dropdown_colmap.addItem("DoLPmap")
        self.dropdown_colmap = dropdown_colmap
        colmap_choice_container.layout().addWidget(dropdown_colmap)

        btn_calculate_colmap = QPushButton("Calculate colourmap")
        btn_calculate_colmap.clicked.connect(self._on_click_btn_calculate_colmap)
                
        
        # histogram plots =====================================================
        graph_container_hist = QWidget()
        graph_container_hist.setMaximumHeight(300)
        graph_container_hist.setLayout(QVBoxLayout())
        
        self.s0_hist_widget = pg.GraphicsLayoutWidget()
        self.s0_hist_widget.setBackground(None)
        graph_container_hist.layout().addWidget(self.s0_hist_widget)
        
        self.dolp_hist_widget = pg.GraphicsLayoutWidget()
        self.dolp_hist_widget.setBackground(None)
        graph_container_hist.layout().addWidget(self.dolp_hist_widget)
        
        
        # layout of the "DoLP limits" gui input ===============================
        dolp_limits_container = QWidget()
        dolp_limits_container.setLayout(QHBoxLayout())
        dolp_limits_container.setMaximumHeight(40)
        # label
        label_dolp_limits = QLabel()
        label_dolp_limits.setText("DoLP min/max limits:")
        label_dolp_limits.setAlignment(Qt.AlignCenter)
        # text edit boxes
        lineedit_dolp_min = QLineEdit()
        lineedit_dolp_max = QLineEdit()
        lineedit_dolp_min.setText("0.0")
        lineedit_dolp_max.setText("1.0")
        ##lineedit_dolp_min.textChanged.connect(self._lineedit_dolp_min_changed)
        ##lineedit_dolp_max.textChanged.connect(self._lineedit_dolp_max_changed)
        self.lineedit_dolp_min = lineedit_dolp_min
        self.lineedit_dolp_max = lineedit_dolp_max
        # combine everything in the container
        dolp_limits_container.layout().addWidget(label_dolp_limits)
        dolp_limits_container.layout().addWidget(lineedit_dolp_min)
        dolp_limits_container.layout().addWidget(lineedit_dolp_max)

        # layout of the "S0 limits" gui input =================================
        s0_limits_container = QWidget()
        s0_limits_container.setLayout(QHBoxLayout())
        s0_limits_container.setMaximumHeight(40)
        # label
        label_s0_limits = QLabel()
        label_s0_limits.setText("S0 min/max limits:")
        label_s0_limits.setAlignment(Qt.AlignCenter)
        # text edit boxes
        lineedit_s0_min = QLineEdit()
        lineedit_s0_max = QLineEdit()
        lineedit_s0_min.setText("0.0")
        lineedit_s0_max.setText("1.0")
        ##lineedit_s0_min.textChanged.connect(self._lineedit_s0_min_changed)
        ##lineedit_s0_max.textChanged.connect(self._lineedit_s0_max_changed)
        self.lineedit_s0_min = lineedit_s0_min
        self.lineedit_s0_max = lineedit_s0_max
        # combine everything in the container
        s0_limits_container.layout().addWidget(label_s0_limits)
        s0_limits_container.layout().addWidget(lineedit_s0_min)
        s0_limits_container.layout().addWidget(lineedit_s0_max)
        
        # layout of "re-render with new limits" gui input =====================
        rerender_colmap_container = QWidget()
        rerender_colmap_container.setLayout(QHBoxLayout())
        rerender_colmap_container.setMaximumHeight(40)
        # label
        label_rerender_colmap = QLabel()
        label_rerender_colmap.setText("Rerender with new min/max limits:")
        label_rerender_colmap.setAlignment(Qt.AlignCenter)
        btn_rerender_colmap = QPushButton("Rerender map")
        btn_rerender_colmap.clicked.connect(self._on_click_rerender_colmap)
        
        # group all colourmap processing gui elements in a box ================
        colmapProcessingGroupBox = QGroupBox("Colourmap rendering")
        basic_processing_box = QVBoxLayout()
        basic_processing_box.addWidget(colmap_choice_container)
        basic_processing_box.addWidget(btn_calculate_colmap)        
        basic_processing_box.addWidget(graph_container_hist)
        basic_processing_box.addWidget(s0_limits_container)
        basic_processing_box.addWidget(dolp_limits_container)
        basic_processing_box.addWidget(btn_rerender_colmap)
        basic_processing_box.addStretch(1)
        colmapProcessingGroupBox.setLayout(basic_processing_box)
        
        
        # =====================================================================
        # Add everything to the layout in order of appearance
        # =====================================================================
        
        self.setLayout(QVBoxLayout())
        #self.layout().setSpacing(0)

        self.layout().addWidget(settingsGroupBox)
        self.layout().addWidget(basicProcessingGroupBox)
        self.layout().addWidget(colmapProcessingGroupBox)
        
        # intialise the histograms
        self.draw_s0_histogram([0,1]) # initialise histogram plot
        self.draw_dolp_histogram([0,1]) # initialise histogram plot

    def check_only_one_layer_is_selected(self):
        numLayersSelected = len(self.viewer.layers.selection)
        if numLayersSelected == 0:
            show_info("Please first select the layer you want to process.")
            return False
        elif numLayersSelected == 1:
            return True
        else:
            show_info(f"Please select only 1 layer. You selected {numLayersSelected} layers.")
            return False
    
    def _on_value_change_bkgnd(self):
        self.offset = float(self.lineedit_bkgnd.text())
    
    def _on_click_quadview(self):
        """" Reorganise the pixels in an unprocessed polarisation camera image
        into a quadview representation. Add as a new layer to the napari viewer.
        Repeat for each layer that was selected. """
        if self.check_only_one_layer_is_selected():
            for layer in self.viewer.layers.selection: # for each selected layer
                pci = PolarisationCameraImage(layer.data,
                                              self.dropdown_method.currentText(),
                                              self.dropdown_unit.currentText(),
                                              self.offset)
                quadview = pci.unprocessed_to_quadview() # calculate quadview
                # add quadview image as new layer to napari viewer
                self.viewer.add_image(quadview, name="quadview_"+layer.name)
                break
            return quadview
    
    def _on_click_channels(self):
        """" Estimate the 4 intensity channels from an unprocessed polarisation
        camera image. Add four new layers (I0, I45, I90 and I135) to the napari
        viewer. Repeat for each layer that was selected. """
        if self.check_only_one_layer_is_selected():
            for layer in self.viewer.layers.selection:
                pci = PolarisationCameraImage(layer.data,
                                              self.dropdown_method.currentText(),
                                              self.dropdown_unit.currentText(),
                                              self.offset)
                I0, I45, I90, I135 = pci.convert_unprocessed()
                I_min = np.min(I0 + I90)/2
                I_max = np.max(I0 + I90)/2
                # add images as new layers to napari viewer
                self.viewer.add_image(I0, contrast_limits=[I_min, I_max], name="I0_"+layer.name)
                self.viewer.add_image(I45, contrast_limits=[I_min, I_max], name="I45_"+layer.name)        
                self.viewer.add_image(I90, contrast_limits=[I_min, I_max], name="I90_"+layer.name)        
                self.viewer.add_image(I135, contrast_limits=[I_min, I_max], name="I135_"+layer.name)
                break
            return I0, I45, I90, I135

    def _on_click_stokes_add_to_viewer(self):
        """" Estimate the Stokes parameter images from an unprocessed polarisation
        camera image. Add three new layers (S0, S1 and S2) to the napari viewer.
        Repeat for each layer that was selected. """
        if self.check_only_one_layer_is_selected():
            for layer in self.viewer.layers.selection:
                
                pci = PolarisationCameraImage(layer.data,
                                              self.dropdown_method.currentText(),
                                              self.dropdown_unit.currentText(),
                                              self.offset)
                pci.subtract_bkgnd()
                I0, I45, I90, I135 = pci.convert_unprocessed()
                
                I0 = I0.astype(np.double)
                I45 = I45.astype(np.double)
                I90 = I90.astype(np.double)
                I135 = I135.astype(np.double)
    
                S0 = (I0 + I45 + I90 + I135)/2
                S1 = I0 - I90
                S2 = I45 - I135
                max_s0 = np.max(S0)
                
                # add a new S0 layer (or replace the data is one was already open)
                if not self.check_if_s0_is_loaded():
                    self.viewer.add_image(S0, contrast_limits=[0, max_s0], name="S0_"+layer.name)    
                else:
                    self.viewer.layers['S0'].data = S0
                # add a new S1 layer (or replace the data is one was already open)
                if not self.check_if_s1_is_loaded():
                    self.viewer.add_image(S1, contrast_limits=[-max_s0, max_s0], name="S1_"+layer.name)    
                else:
                    self.viewer.layers['S1'].data = S1
                # add a new S2 layer (or replace the data is one was already open)
                if not self.check_if_s2_is_loaded():
                    self.viewer.add_image(S2, contrast_limits=[-max_s0, max_s0], name="S2_"+layer.name)    
                else:
                    self.viewer.layers['S2'].data = S2
                    break
                    
            return S0, S1, S2
    
    def calculate_stokes(self):
        """" Estimate the Stokes parameter images from an unprocessed polarisation
        camera image. """
        for layer in self.viewer.layers.selection:
            
            pci = PolarisationCameraImage(layer.data,
                                          self.dropdown_method.currentText(),
                                          self.dropdown_unit.currentText(),
                                          self.offset)
            pci.subtract_bkgnd()
            I0, I45, I90, I135 = pci.convert_unprocessed()
            
            I0 = I0.astype(np.double)
            I45 = I45.astype(np.double)
            I90 = I90.astype(np.double)
            I135 = I135.astype(np.double)

            S0 = (I0 + I45 + I90 + I135)/2
            S1 = I0 - I90
            S2 = I45 - I135
                
        return S0, S1, S2
                     

    def _on_click_aolp(self,add_to_viewer=True):
        """" Calculate the Angle of Linear Polarisation (AoLP) image from the
        S1 and S2 layers that are open in the napari viewer. Add the AoLP image
        to the napari viewer as a new layer."""
        if self.check_only_one_layer_is_selected():
            # calculate Stokes parameters
            if self.checkbox_show_intermediate_results.checkState(): # if checked
                s0, s1, s2 = self._on_click_stokes_add_to_viewer()
            else: # if not checked
                s0, s1, s2 = self.calculate_stokes() # will add as new layers to napari viewer
            
            # calculate AoLP from S1 and S2
            AoLP = (1/2)*np.arctan2(s2,s1)
            # add a new AoLP layer (or replace the data is one was already open)
            if not self.check_if_aolp_is_loaded():
                self.viewer.add_image(AoLP,contrast_limits=[-np.pi/2,np.pi/2])
            else:
                self.viewer.layers['AoLP'].data = AoLP
                
            return AoLP
            
    def _on_click_dolp(self,add_to_viewer=True):
        """" Calculate the Degree of Linear Polarisation (DoLP) image from the
        S0, S1 and S2 layers that are open in the napari viewer. Add the DoLP
        image to the napari viewer as a new layer."""
        # calculate Stokes parameters
        if self.check_only_one_layer_is_selected():
            if self.checkbox_show_intermediate_results.checkState(): # if checked
                s0, s1, s2 = self._on_click_stokes_add_to_viewer()
            else: # if not checked
                s0, s1, s2 = self.calculate_stokes() # will add as new layers to napari viewer
            
            # calculate DoLP from S0, S1 and S2
            DoLP = np.sqrt((s1*s1 + s2*s2)/(s0*s0))
            # add a new DoLP layer (or replace the data is one was already open)
            if not self.check_if_dolp_is_loaded():
                self.viewer.add_image(DoLP,contrast_limits=[0,1])
            else:
                self.viewer.layers['DoLP'].data = DoLP
            
            return DoLP
    
    def _on_click_btn_calculate_colmap(self):
        if self.check_only_one_layer_is_selected():
            if self.dropdown_colmap.currentText() == 'HSVmap':
                self._on_click_hsvmap(None,None)
            elif self.dropdown_colmap.currentText() == 'DoLPmap':
                self._on_click_dolpmap(None,None)
        
    def _on_click_hsvmap(self,lim_s0,lim_dolp):
        # calculate Stokes parameters
        if self.checkbox_show_intermediate_results.checkState(): # if checked
            s0, s1, s2 = self._on_click_stokes_add_to_viewer()
        else: # if not checked
            s0, s1, s2 = self.calculate_stokes() # will add as new layers to napari viewer
        
        if lim_s0 == None: # if no input limits are given
            # update the min/max limits of the s0 box
            s0_min = np.min(s0)
            s0_max = np.max(s0)
            dolp_min = 0.0
            dolp_max = 1.0
        else:
            s0_min = lim_s0[0]
            s0_max = lim_s0[1]
            dolp_min = lim_dolp[0]
            dolp_max = lim_dolp[1]
        
        self.lineedit_s0_min.setText(str(round(s0_min)))
        self.lineedit_s0_max.setText(str(round(s0_max)))
        
        # calculate AoLP and DoLP from Stokes parameters
        AoLP = (1/2)*np.arctan2(s2,s1)
        DoLP = np.sqrt((s1*s1 + s2*s2)/(s0*s0))

        if self.checkbox_show_intermediate_results.checkState(): # if checked
            # add a new AoLP layer (or replace the data is one was already open)
            if not self.check_if_aolp_is_loaded():
                self.viewer.add_image(AoLP,contrast_limits=[-np.pi/2,np.pi/2])
            else:
                self.viewer.layers['AoLP'].data = AoLP
        
            # add a new DoLP layer (or replace the data is one was already open)
            if not self.check_if_dolp_is_loaded():
                self.viewer.add_image(DoLP,contrast_limits=[0,1])
            else:
                self.viewer.layers['DoLP'].data = DoLP
        
        # arrange the hue, saturation and value channels
        h = AoLP # hue
        s = DoLP # saturation
        v = s0 # value
        
        # scale parameters based on limits
        h = (h + np.pi/2)/np.pi; # rescale [-pi/2, pi/2] to [0, 1]
        s = (s - dolp_min)/(dolp_max - dolp_min); # rescale DoLP
        v = (v - s0_min)/(s0_max - s0_min); # rescale DoLP
        # make sure all values fall between 0 and 1
        h[h < 0] = 0
        h[h > 1] = 1
        s[s < 0] = 0
        s[s > 1] = 1
        v[v < 0] = 0
        v[v > 1] = 1
        # convert hsv to rgb
        numDim = len(h.shape) # number of dimensions of the dataset
        hsv = np.stack([h,s,v],numDim) # stack the h, s and v channels along a new dimension
        rgb = hsv_to_rgb(hsv) # convert hsv colourspace to rgb colourspace
        rgb = rgb*255 # scale [0 1] to 8-bit [0 255]
        HSVmap = rgb.astype(np.uint8) # convert to unsigned 8-bit integer values
        HSVmap_R = HSVmap[...,0]
        HSVmap_G = HSVmap[...,1]
        HSVmap_B = HSVmap[...,2]
        
        # add new images for the 3 colour channels: HSVmap_R, HSVmap_G and HSVmap_B
        if 'HSVmap_R' in self.viewer.layers: # if layer already open, replace data
            self.viewer.layers['HSVmap_R'].data = HSVmap_R
        else:
            self.viewer.add_image(HSVmap_R,contrast_limits=[0,255],colormap="red",blending="additive",name="HSVmap_R")

        if 'HSVmap_G' in self.viewer.layers: # if layer already open, replace data
            self.viewer.layers['HSVmap_G'].data = HSVmap_G
        else:
            self.viewer.add_image(HSVmap_G,contrast_limits=[0,255],colormap="green",blending="additive",name="HSVmap_G")

        if 'HSVmap_B' in self.viewer.layers: # if layer already open, replace data
            self.viewer.layers['HSVmap_B'].data = HSVmap_B
        else:
            self.viewer.add_image(HSVmap_B,contrast_limits=[0,255],colormap="blue",blending="additive",name="HSVmap_B")
        
        # link the internal napari contrast limits settings between the layes
        layers_to_link = [self.viewer.layers['HSVmap_R'],self.viewer.layers['HSVmap_G'],self.viewer.layers['HSVmap_B']]
        link_layers(layers_to_link, ('contrast_limits', 'gamma', 'opacity'))
        
        # redraw the histograms
        self.draw_s0_histogram(s0)
        self.draw_dolp_histogram(DoLP)

    def _on_click_dolpmap(self,lim_s0,lim_dolp):
        # calculate Stokes parameters
        if self.checkbox_show_intermediate_results.checkState(): # if checked
            s0, s1, s2 = self._on_click_stokes_add_to_viewer()
        else: # if not checked
            s0, s1, s2 = self.calculate_stokes() # will add as new layers to napari viewer
        
        if lim_s0 == None: # if no input limits are given
            # update the min/max limits of the s0 box
            s0_min = np.min(s0)
            s0_max = np.max(s0)
            dolp_min = 0.0
            dolp_max = 1.0
        else:
            s0_min = lim_s0[0]
            s0_max = lim_s0[1]
            dolp_min = lim_dolp[0]
            dolp_max = lim_dolp[1]
        
        self.lineedit_s0_min.setText(str(round(s0_min)))
        self.lineedit_s0_max.setText(str(round(s0_max)))
        
        # calculate DoLP from S0, S1 and S2
        DoLP = np.sqrt((s1*s1 + s2*s2)/(s0*s0))
        if self.checkbox_show_intermediate_results.checkState(): # if checked
            # add a new DoLP layer (or replace the data is one was already open)
            if not self.check_if_dolp_is_loaded():
                self.viewer.add_image(DoLP,contrast_limits=[0,1])
            else:
                self.viewer.layers['DoLP'].data = DoLP
        
        # scale parameters based on limits
        DoLP = (DoLP - dolp_min)/(dolp_max - dolp_min); # rescale DoLP
        DoLP[DoLP < 0] = 0
        DoLP[DoLP > 1] = 1
        s0 = (s0 - s0_min)/(s0_max - s0_min); # rescale DoLP
        s0[s0 < 0] = 0
        s0[s0 > 1] = 1        
        
        # create rgb stack
        numDim = len(s0.shape) # number of dimensions of the dataset
        r = DoLP
        g = np.ones_like(r)*0.5
        b = 1 - DoLP
        DoLPmap = np.stack([r,g,b],numDim) # stack the h, s and v channels along a new dimension

        # weight all channels by s0
        DoLPmap[...,0] = DoLPmap[...,0]*s0
        DoLPmap[...,1] = DoLPmap[...,1]*s0
        DoLPmap[...,2] = DoLPmap[...,2]*s0
        
        DoLPmap = DoLPmap*255
        DoLPmap = DoLPmap.astype(np.uint8) # convert to unsigned 8-bit integer values
        DoLPmap_R = DoLPmap[...,0]
        DoLPmap_G = DoLPmap[...,1]
        DoLPmap_B = DoLPmap[...,2]
        
        # add new images for the 3 colour channels: HSVmap_R, HSVmap_G and HSVmap_B
        if 'DoLPmap_R' in self.viewer.layers: # if layer already open, replace data
            self.viewer.layers['DoLPmap_R'].data = DoLPmap_R
        else:
            self.viewer.add_image(DoLPmap_R,contrast_limits=[0,255],colormap="red",blending="additive",name="DoLPmap_R")

        if 'DoLPmap_G' in self.viewer.layers: # if layer already open, replace data
            self.viewer.layers['DoLPmap_G'].data = DoLPmap_G
        else:
            self.viewer.add_image(DoLPmap_G,contrast_limits=[0,255],colormap="green",blending="additive",name="DoLPmap_G")

        if 'DoLPmap_B' in self.viewer.layers: # if layer already open, replace data
            self.viewer.layers['DoLPmap_B'].data = DoLPmap_B
        else:
            self.viewer.add_image(DoLPmap_B,contrast_limits=[0,255],colormap="blue",blending="additive",name="DoLPmap_B")
        
        # link the internal napari contrast limits settings between the layes
        layers_to_link = [self.viewer.layers['DoLPmap_R'],self.viewer.layers['DoLPmap_G'],self.viewer.layers['DoLPmap_B']]
        link_layers(layers_to_link, ('contrast_limits', 'gamma', 'opacity'))
        
        # redraw the histograms
        self.draw_s0_histogram(s0)
        self.draw_dolp_histogram(DoLP)
    
    def _on_click_rerender_colmap(self):
        if self.check_only_one_layer_is_selected():
            if self.dropdown_colmap.currentText() == 'HSVmap':
                self._on_click_rerender_hsvmap()
            elif self.dropdown_colmap.currentText() == 'DoLPmap':
                self._on_click_rerender_dolpmap()
        
    def _on_click_rerender_hsvmap(self):
        s0_min = float(self.lineedit_s0_min.text())
        s0_max = float(self.lineedit_s0_max.text())
        dolp_min = float(self.lineedit_dolp_min.text())
        dolp_max = float(self.lineedit_dolp_max.text())
        lim_s0 = (s0_min, s0_max)
        lim_dolp = (dolp_min, dolp_max)
        self._on_click_hsvmap(lim_s0,lim_dolp)
        
    def _on_click_rerender_dolpmap(self):
        s0_min = float(self.lineedit_s0_min.text())
        s0_max = float(self.lineedit_s0_max.text())
        dolp_min = float(self.lineedit_dolp_min.text())
        dolp_max = float(self.lineedit_dolp_max.text())
        lim_s0 = (s0_min, s0_max)
        lim_dolp = (dolp_min, dolp_max)
        self._on_click_dolpmap(lim_s0,lim_dolp)
         
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

    def draw_s0_histogram(self,S0):
        # add a new plot to the s0_hist_widget or empty the old plot
        if not hasattr(self, "plot_s0"):
            self.plot_s0 = self.s0_hist_widget.addPlot()
        else:
            self.plot_s0.clear()

        # draw DoLP histogram
        max_s0 = np.max(S0)
        counts, binedges = np.histogram(S0,bins=100,range=(0.0,max_s0))
        self.plot_s0.plot(binedges[0:-1],counts/np.max(counts),name='S0')
        self.plot_s0.setXRange(0, max_s0, padding=0)
        self.plot_s0.setYRange(0, 1, padding=0)
        self.plot_s0.setLabel('bottom', "S0") # add x-axis label
        
    def draw_dolp_histogram(self,DoLP):
        # add a new plot to the dolp_hist_widget or empty the old plot
        if not hasattr(self, "plot_dolp"):
            self.plot_dolp = self.dolp_hist_widget.addPlot()
        else:
            self.plot_dolp.clear()
        
        # draw DoLP histogram
        counts, binedges = np.histogram(DoLP,bins=100,range=(0.0,1.0))
        self.plot_dolp.plot(binedges[0:-1],counts/np.max(counts),name='DoLP')
        self.plot_dolp.setXRange(0, 1, padding=0)
        self.plot_dolp.setYRange(0, 1, padding=0)
        self.plot_dolp.setLabel('bottom', "DoLP") # add x-axis label



class SetVoxelSize(QWidget):
    
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
        
        self.setLayout(QVBoxLayout())
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
            if len(layer.data.shape) > 2: # if the data has x,y,z dimensions
                layer.scale = [value_z, value_xy, value_xy]
            else: # if the data has no z dimension
                layer.scale = [value_xy, value_xy]