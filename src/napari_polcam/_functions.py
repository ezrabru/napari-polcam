"""
This module contains functions used for polarisation camera image processing.
"""
import numpy as np
from scipy.interpolate import RectBivariateSpline

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
        elif self.method == 'Cubic interpolation':
            I0,I45,I90,I135 = self.interpolate_channels()
        elif self.method == 'Fourier interpolation':
            print('Fourier interpolation is not yet supported.')

        else:
            print('Unexpected value for input variable "method" in method "convert_unprocessed" in class "PolarisationCameraImage".')
        
        return I0, I45, I90, I135
    
    def unprocessed_to_unassigned_channels(self):
        """ Reorganise pixels in unprocessed image to get four unassigned
        intensity channels."""
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
        
        return ch00, ch01, ch10, ch11
    
    def estimate_channels_no_interpolation(self):
        """ Reorganise pixels in unprocessed image to get the 4 intensity
        channels."""
        ch00, ch01, ch10, ch11 = self.unprocessed_to_unassigned_channels()
        # assign channels to correct polariser transmission axis orientation
        I0,I45,I90,I135 = self.assign_channel_to_transmission_axis_orientation(ch00,ch01,ch10,ch11)
        return I0, I45, I90, I135
    
    def interpolate_channels(self):
        
        # get coordinates of the known values in each channel
        ny, nx = self.imgSizeEven
        xx_hr = np.arange(0,nx,1)
        yy_hr = np.arange(0,ny,1)

        # get the known values for each channel (i.e. split channels)
        ch00, ch01, ch10, ch11 = self.unprocessed_to_unassigned_channels()
        
        # low resolution mesh (known mesh) for each channel
        xx00 = np.arange(0,nx,2)
        yy00 = np.arange(0,ny,2)
        
        xx01 = np.arange(1,nx+1,2)
        yy01 = np.arange(0,ny,2)
        
        xx10 = np.arange(0,nx,2)
        yy10 = np.arange(1,ny+1,2)
        
        xx11 = np.arange(1,nx+1,2)
        yy11 = np.arange(1,ny+1,2)
        
        if self.numDim == 2:            
            ch00_interpolated = self.interpolate_frame(xx00, yy00, ch00, xx_hr, yy_hr)
            ch01_interpolated = self.interpolate_frame(xx01, yy01, ch01, xx_hr, yy_hr)
            ch10_interpolated = self.interpolate_frame(xx10, yy10, ch10, xx_hr, yy_hr)
            ch11_interpolated = self.interpolate_frame(xx11, yy11, ch11, xx_hr, yy_hr)
            
        if self.numDim == 3:
            # initialise
            ch00_interpolated = np.zeros_like(self.img)
            ch01_interpolated = np.zeros_like(self.img)
            ch10_interpolated = np.zeros_like(self.img)
            ch11_interpolated = np.zeros_like(self.img)
            numSlices = self.img.shape[0]
            for id_frame in range(numSlices):
                ch00_interpolated[id_frame,:,:] = self.interpolate_frame(xx00, yy00, ch00[id_frame,:,:], xx_hr, yy_hr)
                ch01_interpolated[id_frame,:,:] = self.interpolate_frame(xx01, yy01, ch01[id_frame,:,:], xx_hr, yy_hr)
                ch10_interpolated[id_frame,:,:] = self.interpolate_frame(xx10, yy10, ch10[id_frame,:,:], xx_hr, yy_hr)
                ch11_interpolated[id_frame,:,:] = self.interpolate_frame(xx11, yy11, ch11[id_frame,:,:], xx_hr, yy_hr)
                
        if self.numDim == 4:
            # initialise
            ch00_interpolated = np.zeros_like(self.img)
            ch01_interpolated = np.zeros_like(self.img)
            ch10_interpolated = np.zeros_like(self.img)
            ch11_interpolated = np.zeros_like(self.img)
            numTimepoints = self.img.shape[0]
            numSlices = self.img.shape[1]
            for id_time in range(numTimepoints):
                for id_frame in range(numSlices):
                    ch00_interpolated[id_time,id_frame,:,:] = self.interpolate_frame(xx00, yy00, ch00[id_time,id_frame,:,:], xx_hr, yy_hr)
                    ch01_interpolated[id_time,id_frame,:,:] = self.interpolate_frame(xx01, yy01, ch01[id_time,id_frame,:,:], xx_hr, yy_hr)
                    ch10_interpolated[id_time,id_frame,:,:] = self.interpolate_frame(xx10, yy10, ch10[id_time,id_frame,:,:], xx_hr, yy_hr)
                    ch11_interpolated[id_time,id_frame,:,:] = self.interpolate_frame(xx11, yy11, ch11[id_time,id_frame,:,:], xx_hr, yy_hr)
        
        else:
            print(f"Processing of a {self.numDim}-dimensional dataset is not supported in napari-polcam.")
        
        I0, I45, I90, I135 = self.assign_channel_to_transmission_axis_orientation(ch00_interpolated,ch01_interpolated,ch10_interpolated,ch11_interpolated)
        
        return I0, I45, I90, I135
    
    
    def interpolate_frame(self, xx, yy, ch, xx_hr, yy_hr):
        if self.method == 'Cubic interpolation':
            interp_spline = RectBivariateSpline(xx, yy, ch)
            interpolated_frame = interp_spline(xx_hr, yy_hr)
        #elif self.method == 'Fourier':
            # not yet supported
        else:
            print('Unexpected value for input variable "method" in method "convert_unprocessed" in class "PolarisationCameraImage".')
        
        return interpolated_frame


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
    