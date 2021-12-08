# Packages
import numpy as np 
import sys
sys.path.append('../')
import filters # Filters function used for topology optimization 
import grcwa # Python RCWA Library. See downloading instructions at https://github.com/weiliangjinca/grcwa

args = {}
params = {}

## Define variables 
args['nG'] = 51 # Always check convergence wrt nG 

# Lattice vector (unit length is 1 um)
args['Lx'] = 0.430
args['Ly'] = args['Lx']
args['L1'] = [args['Lx'], 0.]
args['L2'] = [0., args['Ly']]

# Wavelength vector 
wl_vec = np.linspace(0.4,0.8,400)

# Array of z-bound for volume integral (measured depths)
# meas_depth_vec = np.array([100.,200.,300.,400.,500.,600.,700.,800.,900.,1000.])*1e-3 # Coarse sweep
meas_depth_vec = np.linspace(0.01,1.,20) # Thin sweep

# discretization for patterned layers
# Has to be even integers with current C4v implementation
args['Nx'] = 100
args['Ny'] = args['Nx']

## Multilayer setting. Bottom/Top are half spaces, middle is for unpatterned middle layers. 
eps_SiO2 = 2.13353 +1e-5*1j# 
eps_vacuum = 1.

# Geometry
args['thick_max'] = 0.5 # Thickness of silicon top layer 
args['etch_depth'] = 25*1e-3 # To be updated for various hole depths 
thick_SiO2 = 1.0 # Thickness of silica layer 
rad_etch = 0.130 # Radius of hole 

# Angular parameters
theta_max = 17.5*np.pi/180. # Max. angular aperture
theta_res = 20 # Angular resolution 
theta_vec = np.linspace(0, theta_max, theta_res) # Array of angles 

# Two types of polarization
pol_dict = {"s", "p"}

# Loads silicon refractive index (from refractiveindex.info)
eps_Si_data = np.loadtxt(open('permittivity/Si (Silicon) - Palik.txt', "rb"), delimiter=",", skiprows=0)
def eps_Si_func(wl):
    ''' Permittivity of silicon '''
    # wl wavelength in microns
    ind = np.where(np.abs(eps_Si_data[:,0]-wl) == np.min(np.abs(eps_Si_data[:,0]-wl)))       
    return (eps_Si_data[ind,1] + 1j * eps_Si_data[ind,2])[0][0]

def rescaled(dof, dof_min, dof_max):
    ''' Rescales DOF from [0,1] to [dof_min, dof_max'''
    return dof_min + dof*(dof_max-dof_min)

def Veff(freqangs, pol, meas_depth, etch_depth):
    ''' Calculates Veff from experimental structure (Figure 2)
        Given incoming frequency, solid angle (freqangs)
        and polarization (pol) 
        Integrates over measurement depth (meas_depth)
        Circular hole with depth etch_depth
    '''

    # Loads frequency and angles from freqangs
    freq = freqangs[0]
    theta = freqangs[1]
    phi = freqangs[2]

    # Silicon permittivity
    eps_Si = eps_Si_func(1/freq)

    # Creates RCWA Object
    obj = grcwa.obj(args['nG'], args['L1'], args['L2'], freq, theta, phi, verbose = 0)      

    # Creates multi-layer structure 
    obj.Add_LayerUniform(0., eps_vacuum)
    obj.Add_LayerGrid(args['etch_depth'], args['Nx'], args['Ny']) # Layer with circular hole 
    obj.Add_LayerUniform(args['thick_max']-etch_depth, eps_Si)
    obj.Add_LayerUniform(meas_depth, eps_SiO2) # Layer of interest 
    obj.Add_LayerUniform(thick_SiO2 - meas_depth, eps_SiO2) # Rest of silica (not measured here)
    obj.Add_LayerUniform(1000., eps_Si)
    obj.Add_LayerUniform(0., eps_vacuum)

    obj.Init_Setup(Gmethod = 0)
    
    # Updates DOF to make circular hole 
    flattened_dof = 1-filters.dof_to_pillar(rad_etch, args['Nx'], args['Ny'], args['Lx'], args['Ly'], binary_flag = True)
    flattened_dof = np.array(flattened_dof).flatten()*(eps_Si - eps_vacuum) + eps_vacuum # Turns 0/1 DOF into epsilon values        
    obj.GridLayer_geteps(flattened_dof) # Flattens DOFs for RCWA

    # Set incoming polarization
    if pol == "s":
        planewave = {'p_amp': 0, 's_amp': 1, 'p_phase': 0, 's_phase': 0}
    if pol == "p":
        planewave = {'p_amp': 1, 's_amp': 0, 'p_phase': 0, 's_phase': 0}
    
    # Define incoming plane wave 
    obj.MakeExcitationPlanewave(planewave['p_amp'], planewave['p_phase'], planewave['s_amp'], planewave['s_phase'], order = 0)

    SiO2_layer = 3 # index of silica layer of interest 
    dN = 1/args['Nx']/args['Ny'] # Integral discretization
    M0 = grcwa.fft_funs.get_conv(dN, np.ones((args['Nx'], args['Ny'])), obj.G) # Defines integration kernel in Fourier space
    # To be consistent with Poynting vector convention, absorbed power is vol1*omega or vol1/lambda
    # See here: https://github.com/weiliangjinca/grcwa/blob/master/grcwa/rcwa.py
    vol1 = obj.Volume_integral(SiO2_layer, M0, M0, M0, normalize=1) # Calculate volume integral
    
    res = np.abs(vol1) 
    return res

folder_name = 'sweep' # folder to save data in /res/folder_name

# Loop to calculate over all possible parameters. Can be parallelized with multiprocessing.pool.map()
data_mat = np.zeros((len(wl_vec), len(theta_vec), len(theta_vec), len(pol_dict), len(meas_depth_vec))) # Data array 
theta_mat = np.zeros((len(theta_vec), len(theta_vec))) # Theta array (for angular integration)
phi_mat = np.zeros((len(theta_vec), len(theta_vec))) # Phi array (for angular integration)
for (meas_depth_ind, meas_depth) in zip(range(len(meas_depth_vec)), meas_depth_vec):
    args['meas_depth'] = meas_depth
    for (polind, pol) in zip(range(len(pol_dict)), pol_dict):
        for itx in range(len(theta_vec)):
            for ity in range(itx+1):
                # Updates theta and phi
                # Triangular domain sweep within (circular) numerical aperture 
                thetax = theta_vec[itx]
                thetay = theta_vec[ity]
                theta = np.sqrt(thetax**2.+thetay**2.)
                phi = np.arctan(thetay/thetax) if thetax != 0 else np.sign(thetax)*np.pi/2.
                theta_mat[itx, ity] = theta
                phi_mat[itx, ity] = phi
                if theta <= theta_max:
                    # If in angular aperture, calculates RCWA
                    print("Calculating up to depth = {0} nm, pol = {1}, theta = {2}, phi = {3}".format(meas_depth*1e3, pol, theta*180/np.pi, phi*180/np.pi))
                    for (wlind, wl) in zip(range(len(wl_vec)), wl_vec):                        
                        # Calculates RCWA for a given set of parameters
                        freqangs0 = [1/wl, theta, phi]
                        res = Veff(freqangs0, pol, meas_depth, args['etch_depth'])/wl                    
                        data_mat[wlind, itx, ity, polind, meas_depth_ind] = res

# Saves data in big dictionary 
data = {"data_mat":data_mat, "wl_vec": wl_vec, "pol": pol, "theta_vec": theta_vec, "theta_mat": theta_mat, "phi_mat": phi_mat, "meas_depth_vec":meas_depth_vec, "thick_SiO2": thick_SiO2, "rad_etch": rad_etch, "args": args}
# np.save("res/"+folder_name+"/rcwa_expt_analysis_etchdepth_{0}nm".format(args['etch_depth']*1e3), data, allow_pickle = True)