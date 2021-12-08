import sys
import numpy as np 
import matplotlib.pyplot as plt 
import csv
import glob 
import re
import os
from scipy.signal import find_peaks
from scipy.io import loadmat, savemat

class Spectrum():
    '''
    Creates spectrum object which contains 
        - wavelength array
        - list of intensity array
        - dictionary with information on the measurement (current, beam voltage, integration time, etc.)

        ## System of Units ##
        Current is in uA 
        Beam voltage is in keV
        Integration time is in s 
    '''
    def __init__(self, wl, lint, name, dict = {}):
        self.wl = wl
        self.lint = lint
        self.err = [0.]
        self.name = name
        if dict == {}:
            self.dict = {"curr": 1, "volt": 40, "intime": 1}
        else:
            self.dict = dict

    def add_info(self, key, dict_entry):
        self.dict["key"] = dict_entry

    def add_data(self, new_spectrum):
        if len(self.wl) == 0:
            self.wl = new_spectrum.wl
            self.lint = new_spectrum.lint
        else:
            if not (self.wl == new_spectrum.wl).all():
                raise ValueError("Wavelength domain of measurements do not match. Check consistency of the data.")
            self.lint = np.vstack((self.lint, new_spectrum.lint))

    def remove_pixel(self, pixel_ind, verbose = True):
        # Removes pixel in data given index 
        # Use when single first pixel is dead     
        if verbose:
            print("Removing manually pixel {0} from data with wl = {1} signal = {2}".format(pixel_ind, self.wl[pixel_ind], self.lint[pixel_ind]))
        self.wl = np.delete(self.wl, pixel_ind)
        self.lint = np.delete(self.lint, pixel_ind)
        if len(self.err)>1:
            self.err = np.delete(self.err, pixel_ind)

    def remove_background(self, background_spectrum):
        # Assumes background is already averaged (signal = 1d array) 
        if len(self.wl) == 0:
            self.wl = background_spectrum.wl
            self.lint = -background_spectrum.lint
        else:
            if not (self.wl == background_spectrum.wl).all():
                raise ValueError("Wavelength domain of background and measurements do not match. Check consistency of the data.")
            if len(np.shape(self.lint)) > 1:
                for ind in range(np.shape(self.lint)[1]):
                    self.lint[ind, :] -= background_spectrum.lint
            else:
                self.lint -= background_spectrum.lint

    def remove_cosmic(self, threshold = 2.):
        # Removes spikes from cosmic radiation from the spectrum
        peaks, properties = find_peaks(self.lint, threshold = threshold, height=(None,None), prominence=(15.,None), width=(None,3))       
        # peaks, properties = find_peaks(self.lint, threshold = None, height=(250,None), prominence=(None,None), width=(None,None))       
#         print("Peak properties: nb of peaks = {} heights = {} prominences = {} widths = {}".format(len(properties["peak_heights"]), properties["peak_heights"], properties["prominences"], properties["widths"]))
        if len(properties["peak_heights"]) > 1:
            self.lint[peaks-1] = self.lint[peaks-2]
            self.lint[peaks] = self.lint[peaks+2]
            self.lint[peaks+1] = self.lint[peaks+2]
            self.err[peaks-1] = self.err[peaks-2]
            self.err[peaks] = self.err[peaks+2]
            self.err[peaks+1] = self.err[peaks+2]

    def add_error(self, error):
        # Adds constant error (if scalar) or errorbar array (otherwise) to data
        if len(error) == 1:
            # Assume standard deviation is wavelength-invariant
            self.err = np.ones_like(self.lint)*error
        else:
            if len(error) != len(self.lint):
                raise ValueError("Errors and intensities should have the same size.")
            else:
                self.err = error

    def process_spectra_list(self):
        # When lint has more than a single element, 
        # computes mean and standard deviation of data
        if len(np.shape(self.lint)) == 1:
            print("Single measurement, no error can be computed from data. Add manually with self.add_error.")
        else:
            self.err = np.std(self.lint, axis = 0)
            self.lint = np.mean(self.lint, axis = 0)

    def power_calibrate(self, loss_function):
        # Normalizes signal and error to a known loss function
        self.lint = np.divide(self.lint, loss_function)
        self.err = np.divide(self.err, loss_function)

    def plot_spectra(self, fig, normalization = "current", color = ""):
        # fig must be matplotlib figure object 
        if len(self.err) > 0 or self.err[0] != 0:
            if normalization == "current":
                if color != "":
                    p = plt.plot(self.wl, self.lint/self.dict["curr"], label = self.dict["caption"], color = color)
                else:
                    p = plt.plot(self.wl, self.lint/self.dict["curr"], label = self.dict["caption"])
                prev_color = p[0].get_color()
                p = plt.fill_between(self.wl, (self.lint+self.err)/self.dict["curr"], (self.lint-self.err)/self.dict["curr"], alpha = 0.3, color = prev_color)
            if normalization == "none":
                if color != "":
                    p = plt.plot(self.wl, self.lint/self.dict["curr"], label = self.dict["caption"], color = color)
                else:
                    p = plt.plot(self.wl, self.lint/self.dict["curr"], label = self.dict["caption"])
                prev_color = p[0].get_color()
                p = plt.fill_between(self.wl, self.lint+self.err, self.lint-self.err, alpha = 0.3, color = prev_color)                
        else:
            p = plt.plot(self.wl, self.lint)
        return fig 

    def save_spectra(self, file_name = "", extension = ".npy"):
        data = {"wl": self.wl, "lint": self.lint, "err": self.err, "name": self.name, "dict": self.dict}
        if extension == ".npy":
            if file_name == "":
                np.save(self.name + ".npy", data)
            else:
                np.save(file_name, data, allow_pickle=True)
        if extension == ".mat":
            if file_name == "":
                savemat(self.name + ".mat", data)
            else:
                savemat(file_name, data)

def read_data_file(data_file, extension = 'csv'):    
    wl = []
    intensity = []
    if extension == 'csv':
        with open(data_file) as csvDataFile:
            csvReader = csv.reader(csvDataFile)
            next(csvReader)
            for row in csvReader:
                wl.append(float(row[0]))
                intensity.append(float(row[1]))
    if extension == 'txt':
        # Specific to txt file
        # Format from old software (Windows XP)
        data_mat = np.loadtxt(data_file) 
        wl = data_mat[:, 0]
        intensity = data_mat[:, 2]
    return Spectrum(np.array(wl), np.array(intensity), "")

def find_all_spectra(folder, caption_flag = 'CLPhC/', caption_length = 6, extension = 'csv', sorted_flag = False, split_flag = '-Frame-[0-9]+'):
    # Returns a list of all spectra in a folder 
    # Finds spectra that are Frames of a longer integration time
    # and puts them together in a single Spectrum    
    if not(sorted_flag):
        all_csv_files = glob.glob(os.getcwd() + folder + '/*.{}'.format(extension))
    if sorted_flag:
        all_csv_files = sorted(glob.glob(os.getcwd() + folder + '/*.{}'.format(extension)))
    all_unique_meas = [re.split(split_flag, csv_file)[0] for csv_file in all_csv_files]
    all_unique_meas = np.unique(all_unique_meas)
    
    spectra = []
    for meas in all_unique_meas:
        count_spectra = 0 
        print(caption_flag, meas)
        cap = re.split(caption_flag, meas)[1]
        cap = cap[:caption_length]
        spec_temp = Spectrum([], [], meas, {}) 
        spec_temp.dict["caption"]  = cap        
        # Find all csv files corresponding to this measurement 
        for csv_file in all_csv_files:
            if (meas in csv_file) or meas == csv_file:
                count_spectra += 1
                spec_csv = read_data_file(csv_file, extension)
                wl_file = spec_csv.wl
                int_file = spec_csv.lint
                spec_temp.add_data(spec_csv)
        spec_temp.dict["nframes"] = count_spectra
        spectra.append(spec_temp)
        print("Measurement {}, {} integrated files".format(meas, count_spectra))
    return spectra

def find_all_current_files_match(folder, list_of_spectra, caption_flag = "without_polarizer", extension = "mat"):
    '''
    Reads all Matlab .mat files in folder and matches them to existing spectra in list 
    caption_flag allows to capture only the relevant part of the file's name
    Should be run after reading all spectra and compiling them into single files (accumulations)
    '''
    all_mat_files = glob.glob(os.getcwd() + folder + '/*.{}'.format(extension))
    for current_data in all_mat_files:
        # Two types of caption writing depending on file type (1st works for Lightfield, 2nd for Windows XP Software)
        cap1 = re.split(caption_flag, current_data)[1][1:-4]+" "
        # Match to spectrum in list 
        for spec in list_of_spectra:
            if re.search(cap1, spec.name):
                # Open file and read current data                 
                loaded_current_data = loadmat(current_data)
                # Add to matching spectrum 
                spec.dict["curr"] = loaded_current_data["meacurr"][0][0]
                spec.dict["curr_std"] = loaded_current_data["stdcurr"][0][0]              

    return list_of_spectra