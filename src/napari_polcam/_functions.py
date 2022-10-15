"""
This module contains functions used for polarisation camera image processing.
"""
import numpy as np


class PolarisationCameraImage():
    
    def __init__(self, img, method, polariser_unit):
        super().__init__()
        
        self.method = method
        self.polariser_unit = polariser_unit
        
        self.imgSizeEven = (img.shape[0], img.shape[1])
        
        self.mask0
        self.mask45
        self.mask90
        self.mask135
    
    def convert_unprocessed(self):
        if self.method == 'linear interpolation':  
            print('Linear interpolation')
        elif self.method == 'cubic interpolation':
            print('Cubic interpolation')
        elif self.method == 'cubic spline interpolation':
            print('Cubic spline interpolation')
        elif self.method == 'fourier':
            print('Fourier interpolation')
        elif self.method == 'none':
            print('No interpolation')
        else:
            print('Unexpected value for input variable "method" in method "convert_unprocessed" in class "PolarisationCameraImage".')
    
    def estimate_stokes_no_interpolation(self):
        print('No interpolation')
        
    def estimate_stokes_linear_interpolation(self):
        print('Linear interpolation')
        
    def estimate_stokes_cubic_interpolation(self):
        print('Cubic interpolation')
        
    def estimate_stokes_cubic_spline_interpolation(self):
        print('Cubic spline interpolation')
        
    def estimate_stokes_fourierinterpolation(self):
        print('Fourier interpolation')
    
    def aolp_from_stokes(S1,S2):
        return (1/2)*np.arctan(S2/S1)
    
    def dolp_from_stokes(S0,S1,S2):
        return np.sqrt((S1^2 + S2^2)/(S0^2))
    
    def intensities_from_stokes(S0,S1,S2):
        I0   = (S0 + S1)/2;
        I90  = (S0 - S1)/2;
        I45  = (S0 + S2)/2;
        I135 = (S0 - S2)/2;
        return I0, I45, I90, I135
    
    def stokes_from_intensities(I0,I45,I90,I135):
        S0 = (I0 + I45 + I90 + I135)/2
        S1 = I0 - I90
        S2 = I45 - I135
        return S0, S1, S2
    
    def get_masks_polarizer_array(self):
        """ Get binary masks for the differently
        %oriented polarisers in the image plane.
        """
        # generate the masks for the correct number of pixels or slightly more
        numRows, numCols = self.imgSizeEven
        numRepsMosaic = np.ceil(np.max(numRows,numCols)/2);
        [mask00,mask01,mask10,mask11] = self.get_mosaic_masks(numRepsMosaic);
                    
        # remove any excess rows or columns
        mask00 = mask00[0:numRows+1][0:numCols+1]
        mask01 = mask01[0:numRows+1][0:numCols+1]
        mask10 = mask10[0:numRows+1][0:numCols+1]
        mask11 = mask11[0:numRows+1][0:numCols+1]
        
        # assign masks to right pixels
        if self.polariser_unit[0][0] == 0: # top left pixel of image has a polarizer with transmission axis at 0 degrees
            self.mask0  = mask00
            self.mask90 = mask11
            if self.polariser_unit[0][1] == 45:
                self.mask45 = mask01
                self.mask135 = mask10
            else:
                self.mask45 = mask10
                self.mask135 = mask01

        elif self.polariser_unit[0][0] == -45: # top left pixel of image has a polarizer with transmission axis at -45 degrees
            self.maskn45 = mask00
            self.maskp45 = mask11
            if self.polariser_unit[0][1] == 0:
                self.mask0  = mask01
                self.mask90 = mask10
            else:
                self.mask0  = mask10
                self.mask90 = mask01

        elif self.polariser_unit[0][0] == 45: # top left pixel of image has a polarizer with transmission axis at 45 degrees
            self.maskp45 = mask00
            self.maskn45 = mask11
            if self.polariser_unit[0][1] == 0:
                self.mask0  = mask01
                self.mask90 = mask10
            else:
                self.mask0  = mask10
                self.mask90 = mask01

        elif self.polariser_unit[0][0] == 90: # top left pixel of image has a polarizer with transmission axis at 90 degrees
            self.mask90 = mask00
            self.mask0  = mask11
            if self.polariser_unit[0][1] == 45:
                self.mask45 = mask01
                self.mask135 = mask10
            else:
                self.mask45 = mask10
                self.mask135 = mask01
        else:       
            print('Unexpected value for input variable "polariser_unit" in method "get_masks_polarizer_array" in class "PolarisationCameraImage".')
    
    
    def get_mosaic_masks(numRepsMosaic):
        mask00 = np.array([[1, 0], [0, 0]])
        mask01 = np.array([[0, 1], [0, 0]])
        mask10 = np.array([[0, 0], [1, 0]])
        mask11 = np.array([[0, 0], [0, 1]])
        mask00 = np.tile(mask00, (numRepsMosaic,numRepsMosaic))
        mask01 = np.tile(mask01, (numRepsMosaic,numRepsMosaic))
        mask10 = np.tile(mask10, (numRepsMosaic,numRepsMosaic))
        mask11 = np.tile(mask11, (numRepsMosaic,numRepsMosaic))
        return mask00, mask01, mask10, mask11
    
    
    
    
    
    
    
    
    
    