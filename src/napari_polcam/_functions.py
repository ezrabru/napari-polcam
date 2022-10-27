"""
This module contains functions used for polarisation camera image processing.
"""
import numpy as np
from scipy.interpolate import interp2d, RectBivariateSpline


class PolarisationCameraImage():
    
    def __init__(self, img, method, polariser_unit, offset):
        super().__init__()
        
        self.img = img - offset
        self.method = method
        self.polariser_unit = polariser_unit
        self.offset = offset
                
        self.numDim = np.ndim(self.img) # number of dimensions of the dataset
        
        nrows = img.shape[-1]
        ncols = img.shape[-2]
        if np.mod(nrows,2):
            nrows = nrows - 1
        if np.mod(ncols,2):
            ncols = ncols - 1
        self.imgSizeEven = (nrows,ncols)
        self.make_dimensions_even() # make image dimensions in x and y even
    
    def subtract_bkgnd(self):
        bkgnd_corrected_img = self.img - self.offset
        bkgnd_corrected_img[bkgnd_corrected_img < 0] = 0
        self.img = bkgnd_corrected_img
    
    def unprocessed_to_quadview(self):

        if self.numDim == 2:
            # demosaick the four channels
            ch00 = self.img[::2, ::2]
            ch01 = self.img[::2, 1::2]
            ch10 = self.img[1::2, ::2]
            ch11 = self.img[1::2, 1::2]
            # tile the channels together
            quadview_top = np.concatenate((ch00,ch01),axis=1)
            quadview_btm = np.concatenate((ch10,ch11),axis=1)
            quadview = np.concatenate((quadview_top,quadview_btm),axis=0)

        elif self.numDim == 3:
            # demosaick the four channels
            ch00 = self.img[:, ::2, ::2]
            ch01 = self.img[:, ::2, 1::2]
            ch10 = self.img[:, 1::2, ::2]
            ch11 = self.img[:, 1::2, 1::2]
            # tile the channels together
            quadview_top = np.concatenate((ch00,ch01),axis=2)
            quadview_btm = np.concatenate((ch10,ch11),axis=2)
            quadview = np.concatenate((quadview_top,quadview_btm),axis=1)
        
        elif self.numDim == 4:
            # demosaick the four channels
            ch00 = self.img[:, :, ::2, ::2]
            ch01 = self.img[:, :, ::2, 1::2]
            ch10 = self.img[:, :, 1::2, ::2]
            ch11 = self.img[:, :, 1::2, 1::2]
            # tile the channels together
            quadview_top = np.concatenate((ch00,ch01),axis=3)
            quadview_btm = np.concatenate((ch10,ch11),axis=3)
            quadview = np.concatenate((quadview_top,quadview_btm),axis=2)
        
        else:
            quadview = None
        
        return quadview
    
    def convert_unprocessed(self):
        if self.method == 'None':
            I0,I45,I90,I135 = self.estimate_channels_no_interpolation()
        elif self.method == 'Linear interpolation':  
            I0,I45,I90,I135 = self.estimate_channels_linear_interpolation()
        elif self.method == 'Cubic interpolation':
            I0,I45,I90,I135 = self.estimate_channels_cubic_interpolation()
        elif self.method == 'Bivariate spline interpolation':
            I0,I45,I90,I135 = self.estimate_channels_spline_interpolation()
        elif self.method == 'Fourier':
            I0,I45,I90,I135 = self.estimate_channels_fourier_interpolation()
        else:
            print('Unexpected value for input variable "method" in method "convert_unprocessed" in class "PolarisationCameraImage".')
        
        return I0, I45, I90, I135
    
    def estimate_channels_no_interpolation(self):
        """ Reorganise pixels in unprocessed image to get the 4 intensity
        channels."""
        if self.numDim == 2:
            # demosaick the four channels
            ch00 = self.img[::2, ::2]
            ch01 = self.img[::2, 1::2]
            ch10 = self.img[1::2, ::2]
            ch11 = self.img[1::2, 1::2]

        elif self.numDim == 3:
            # demosaick the four channels
            ch00 = self.img[:, ::2, ::2]
            ch01 = self.img[:, ::2, 1::2]
            ch10 = self.img[:, 1::2, ::2]
            ch11 = self.img[:, 1::2, 1::2]
        
        elif self.numDim == 4:
            # demosaick the four channels
            ch00 = self.img[:, :, ::2, ::2]
            ch01 = self.img[:, :, ::2, 1::2]
            ch10 = self.img[:, :, 1::2, ::2]
            ch11 = self.img[:, :, 1::2, 1::2]
        
        else:
            print(f"Processing of a {self.numDim}-dimensional dataset is not supported in napari-polcam.")
        
        # assign channels to correct polariser transmission axis orientation
        I0,I45,I90,I135 = self.assign_channel_to_transmission_axis_orientation(ch00,ch01,ch10,ch11)
        return I0, I45, I90, I135
        
    def estimate_channels_linear_interpolation(self):
        I0 = self.interpolate_channel(self.img,self.mask0,'linear')
        I45 = self.interpolate_channel(self.img,self.mask45,'linear')
        I90 = self.interpolate_channel(self.img,self.mask90,'linear')
        I135 = self.interpolate_channel(self.img,self.mask135,'linear')
        return I0,I45,I90,I135
        
    def estimate_channels_cubic_interpolation(self):
        I0 = self.interpolate_channel(self.img,self.mask0,'cubic')
        I45 = self.interpolate_channel(self.img,self.mask45,'cubic')
        I90 = self.interpolate_channel(self.img,self.mask90,'cubic')
        I135 = self.interpolate_channel(self.img,self.mask135,'cubic')
        return I0,I45,I90,I135
        
    def estimate_channels_spline_interpolation(self):
        RectBivariateSpline
        I0 = self.interpolate_channel_spline(self.img,self.mask0)
        I45 = self.interpolate_channel_spline(self.img,self.mask45)
        I90 = self.interpolate_channel_spline(self.img,self.mask90)
        I135 = self.interpolate_channel_spline(self.img,self.mask135)
        return I0,I45,I90,I135
        
    def estimate_channels_fourier_interpolation(self):
        print('Fourier interpolation')
    
    def interpolate_channel(self,img,mask,interp2d_method):
        nx, ny = self.imgSizeEven
        
        if mask[0][0]:
            xx,yy = np.meshgrid(np.arange(1,nx+1,2), np.arange(1,ny+1,2)) 
        elif mask[0][1]:
            xx,yy = np.meshgrid(np.arange(2,nx+2,2), np.arange(1,ny+1,2)) 
        elif mask[1][0]:
            xx,yy = np.meshgrid(np.arange(1,nx+1,2), np.arange(1,ny+2,2)) 
        elif mask[1][1]:
            xx,yy = np.meshgrid(np.arange(2,nx+2,2), np.arange(1,ny+2,2)) 
        else:
            print("Unexpected mask in function 'interpolate_channel' in class 'PolarisationCameraImage'.")

        z = img*mask # reshape image first
        f = interp2d(xx, yy, z, kind='cubic')
        xnew, ynew = np.meshgrid(np.arange(1,nx+1,1), np.arange(1,ny+1,1))
        return f(xnew, ynew)
    
    def interpolate_channel_spline(self,img,mask):
        nx, ny = self.imgSizeEven
        
        if mask[0][0]:
            xx,yy = np.meshgrid(np.arange(1,nx+1,2), np.arange(1,ny+1,2)) 
        elif mask[0][1]:
            xx,yy = np.meshgrid(np.arange(2,nx+2,2), np.arange(1,ny+1,2)) 
        elif mask[1][0]:
            xx,yy = np.meshgrid(np.arange(1,nx+1,2), np.arange(1,ny+2,2)) 
        elif mask[1][1]:
            xx,yy = np.meshgrid(np.arange(2,nx+2,2), np.arange(1,ny+2,2)) 
        else:
            print("Unexpected mask in function 'interpolate_channel_spline' in class 'PolarisationCameraImage'.")

        z = img*mask # reshape image first
        f = RectBivariateSpline(xx, yy, z)
        xnew, ynew = np.meshgrid(np.arange(1,nx+1,1), np.arange(1,ny+1,1))
        return f(xnew, ynew)
    
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
    
    def assign_channel_to_transmission_axis_orientation(self,ch00,ch01,ch10,ch11):
        if self.polariser_unit == "[-45 0; 90 45]":
            I0 = ch01
            I45 = ch11
            I90 = ch10
            I135 = ch00
        elif self.polariser_unit == "[0 -45; 45 90]":
            I0 = ch00
            I45 = ch10
            I90 = ch11
            I135 = ch01
        elif self.polariser_unit == "[90 45; -45 0]":
            I0 = ch11
            I45 = ch01
            I90 = ch00
            I135 = ch10
        elif self.polariser_unit == "[45 90; 0 -45]":
            I0 = ch10
            I45 = ch00
            I90 = ch01
            I135 = ch11
        else:
            print(f"The value {self.polariser_unit} is not supported.")
                
        return I0, I45, I90, I135
    
    def make_dimensions_even(self):
        nrows, ncols = self.imgSizeEven
        if self.numDim == 2:
            self.img = self.img[:ncols, :nrows]
        elif self.numDim == 3:
            self.img = self.img[:, :ncols, :nrows]
        elif self.numDim == 4:
            self.img = self.img[:, :, :ncols, :nrows]        
        else:
            print(f"Processing of a {self.numDim}-dimensional dataset is not supported in napari-polcam.")

    
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
    
    def aolp_from_stokes(S1,S2):
        return (1/2)*np.arctan(S2/S1)
    
    def dolp_from_stokes(S0,S1,S2):
        return np.sqrt((S1**2 + S2**2)/(S0**2))
    
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
    