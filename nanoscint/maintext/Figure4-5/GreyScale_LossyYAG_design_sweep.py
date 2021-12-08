# Packages
import numpy as np
import sys
sys.path.append('../')
import permittivities
from utils import rescaled
import filters # Filters function used for topology optimization 
import grcwa # Python RCWA Library. See downloading instructions at https://github.com/weiliangjinca/grcwa

def runsweep(period, loss):
    ''' Runs Veff sweep for a given period of the photonic crystal and loss of YAG:Ce
        All length units in microns
        Loss in units of permittivity (imaginary)
    '''
    args = {}

    ## Define variables 
    args['nG'] = 51 # Always check convergence wrt nG

    # Define square attice vector (unit length is 1 um)
    args['Lx'] = period
    args['Ly'] = args['Lx']
    args['L1'] = [args['Lx'], 0.]
    args['L2'] = [0., args['Ly']]

    # Calculated wavelengths
    wl_vec = np.linspace(0.540, 0.560, 500)
            
    # Discretization for patterned layers
    # Has to be even integers with current C4v implementation
    args['Nx'] = 200
    args['Ny'] = args['Nx']

    # Geometry
    args['rad_max'] = args['Lx']/2.    

    ''' Pick experimental design
        thick_max is thickness of YAG:Ce substrate        
        etch_depth is depth of hole (units of thick_max)
        etch_rad is radius of hole, etch_rad = 0 makes reference simulation (units of rad_max)
        "Greyscale" simulation: etch_rad not needed (size of hole set by sin^2 profile)
        Binary simulation: etch_rad defines size of hole 
        Comment appropriate line in Veff function to switch between binary/greyscale
    '''
    # # Dose A2, 50\mu m (Figure 5)
    args['thick_max'] = 50.0
    etch_depth_vec = np.array([0.0340])/args['thick_max'] # Measured depth 
    # etch_depth_vec = np.array([0.034,0.024,0.044])/args['thick_max'] # Measured depth +/- uncertainty
    etch_rad_vec = np.array([0., 0.200/2/args['rad_max']])             

    # Dose A1, 20\mu m (Figure 4)
    # args['thick_max'] = 20.0
    # etch_depth_vec = np.array([0.050/args['thick_max']])
    # etch_depth_vec = np.array([0.04,0.050,0.060])/args['thick_max']  
    # etch_rad_vec = np.array([0., 0.200/2/args['rad_max']])                 

    # Angular sweep 
    theta_max = 17.5*np.pi/180. # half-angle numerical aperture 
    theta_res = 8 # resolution along one axis 
    theta_vec = np.linspace(0, theta_max, theta_res)

    pol_dict = {"s", "p"}
    eps_vacuum = 1.0

    def Veff(freqangs, dof, pol, ref_flag = False, method_flag = "Vol"):
        ''' Calculates Veff from experimental structure (Figure 4-5)
        Given incoming frequency, solid angle (freqangs)
        and polarization (pol) 
        Integrates over measurement depth (meas_depth)
        Circular hole with depth etch_depth
        if ref_flag: reference simulation with no etched hole 
        Available methods: 
            - "Vol" (Calculates Volume integrale)
            - "Abs" (Calculates absorption as 1 - R - T): only works if scintillator is the only absorbing material
        '''
        # dof[0] = etch depth
        # dof[1] = hole radius 
        freq = freqangs[0]
        theta = freqangs[1]
        phi = freqangs[2]

        eps_scint = np.real(permittivities.eps_YAG(1/freq)) + loss*1j

        # Define etch depth 
        etch_depth = rescaled(dof[0], 0., args['thick_max'])

        # Define RCWA object 
        obj = grcwa.obj(args['nG'], args['L1'], args['L2'], freq, theta, phi, verbose = 0)      
        obj.Add_LayerUniform(0., eps_vacuum)
        if not(ref_flag):
            # Simulation with etched hole 
            ## Etch first, uniform then
            loi_pat = 1 # layer of interest (pat)
            loi_unpat = 2 # layer of interest (unpat)
            obj.Add_LayerGrid(etch_depth, args['Nx'], args['Ny'])
            obj.Add_LayerUniform(args['thick_max']-etch_depth, eps_scint)        
        if ref_flag:
            # Reference simulation: bulk 
            obj.Add_LayerUniform(args['thick_max'], eps_scint)
        obj.Add_LayerUniform(0., eps_vacuum)
        obj.Init_Setup(Gmethod = 0)
        
        if not(ref_flag):
            ''' This is where one can switch between binary / greyscale DOFs '''
            # non_flattened_dof = 1. - filters.dof_to_pillar(hole_radius, args['Nx'], args['Ny'], args['Lx'], args['Ly'], binary_flag=False) # Binary
            non_flattened_dof = 1. - filters.sin_profile(args['Nx'], args['Ny'], args['Lx'], args['Ly']) # Greyscale 
            # Flatten DOF for RCWA
            flattened_dof = non_flattened_dof.flatten()
            flattened_dof = np.array(flattened_dof).flatten()*(eps_scint - eps_vacuum) + eps_vacuum
            obj.GridLayer_geteps(flattened_dof)

        # Defines incoming polarization and plane wave object
        if pol == "s":
            planewave = {'p_amp': 0, 's_amp': 1, 'p_phase': 0, 's_phase': 0}
        if pol == "p":
            planewave = {'p_amp': 1, 's_amp': 0, 'p_phase': 0, 's_phase': 0}
        obj.MakeExcitationPlanewave(planewave['p_amp'], planewave['p_phase'], planewave['s_amp'], planewave['s_phase'], order = 0)
        
        # Calculates reflection-transmission from multilayer
        R, T = obj.RT_Solve(normalize = 1)
        if method_flag == "Vol":
            # Volume integrale method
            dN = 1/args['Nx']/args['Ny']
            if not(ref_flag):
                M1 = grcwa.fft_funs.get_conv(dN, non_flattened_dof, obj.G)
                vol1 = obj.Volume_integral(loi_pat, M1, M1, M1, normalize=1)
                
                M2 = grcwa.fft_funs.get_conv(dN, np.ones((args['Nx'], args['Ny'])), obj.G)
                vol2 = obj.Volume_integral(loi_unpat, M2, M2, M2, normalize=1)

                res = np.abs(vol1 + vol2) 
            if ref_flag:
                M1 = grcwa.fft_funs.get_conv(dN, np.ones((args['Nx'], args['Ny'])), obj.G)
                vol1 = obj.Volume_integral(1, M1, M1, M1, normalize=1)

                res = np.abs(vol1)

        if method_flag == "Abs":        
            # Absorption method
            res = 1.0-R-T

        return (res, R, T)

    folder_name = 'new_sweep' # Data is saved in /res/folder_name

    ''' Loop to calculate over all parameters. Can be parallelized with multiprocessing.pool.map() '''
    # Updates etch depth (normalized DOF)
    for (etch_depth_ind, etch_depth_dof) in zip(range(len(etch_depth_vec)), etch_depth_vec):
        etch_depth_rescaled = rescaled(etch_depth_dof, 0., args['thick_max'])*1e3
        # Updates etch radius (normalized DOF)
        for (etch_rad_ind, etch_rad_dof) in zip(range(len(etch_rad_vec)), etch_rad_vec):
            etch_rad_rescaled = rescaled(etch_rad_dof, 0., args['rad_max'])*1e3
            dof = np.array([etch_depth_dof, etch_rad_dof])
            # Initializes data and angular sweep data 
            data_mat = np.zeros((len(wl_vec), len(theta_vec), len(theta_vec), len(pol_dict)))
            R_mat = np.zeros((len(wl_vec), len(theta_vec), len(theta_vec), len(pol_dict)))
            T_mat = np.zeros((len(wl_vec), len(theta_vec), len(theta_vec), len(pol_dict)))
            theta_mat = np.zeros((len(theta_vec), len(theta_vec)))
            phi_mat = np.zeros((len(theta_vec), len(theta_vec)))
            # Updates incident polarization
            for (polind, pol) in zip(range(len(pol_dict)), pol_dict):
                count_angles = 0 
                # Updates incident angles 
                # Triangular domain sweep within (circular) numerical aperture
                for itx in range(len(theta_vec)):
                    for ity in range(itx+1):
                        thetax = theta_vec[itx]
                        thetay = theta_vec[ity]
                        theta = np.sqrt(thetax**2.+thetay**2.)
                        phi = np.arctan(thetay/thetax) if thetax != 0 else np.sign(thetax)*np.pi/2.
                        theta_mat[itx, ity] = theta
                        phi_mat[itx, ity] = phi
                        if theta <= theta_max:
                            print("Cal. L = {6}nm, etch depth = {3}nm/{5}nm, rad = {4}nm, pol = {0}, theta = {1}, phi = {2}, losses = {7}".format(pol, theta*180/np.pi, phi*180/np.pi, etch_depth_rescaled, etch_rad_rescaled, args["thick_max"]*1e3, args['Lx']*1e3, loss))
                            count_angles += 1
                            for (wlind, wl) in zip(range(len(wl_vec)), wl_vec):                        
                                freqangs0 = [1/wl, theta, phi]
                                # print("Calculating wl {0} nm".format(wl*1e3))                
                                if etch_rad_rescaled == 0:
                                    # res, R, T = Veff(freqangs0, dof, pol, ref_flag=True, method_flag="Vol")
                                    res, R, T =  Veff(freqangs0, dof, pol, ref_flag=True, method_flag="Abs")
                                else:
                                    # res, R, T = Veff(freqangs0, dof, pol, ref_flag=False, method_flag="Vol")
                                    res, R, T =  Veff(freqangs0, dof, pol, ref_flag=False, method_flag="Abs")
                                data_mat[wlind, itx, ity, polind] = res
                                R_mat[wlind, itx, ity, polind] = R
                                T_mat[wlind, itx, ity, polind] = T

            # Saves data for every (etch_depth, etch_radius) 
            print("Total number of angles {0}".format(count_angles))
            data = {"data_mat":data_mat, "T_mat":T_mat, "R_mat":R_mat, "wl_vec": wl_vec, "pol": pol, "theta_vec": theta_vec, 
            "theta_mat": theta_mat, "phi_mat": phi_mat, "theta_max": theta_max, "etch_rad": etch_rad_rescaled, 
            "etch_depth": etch_depth_rescaled, "etch_depth_vec": etch_depth_vec, "etch_rad_vec": etch_rad_vec, 
            "args": args, "loss": loss}
            # Save data 
            np.save("res/"+folder_name+"/rcwa_YAG_{5}um_design_L_{3}nm_etch_depth_{0}nm_radius{1}nm_nG_{2}_loss_{4}".format(etch_depth_rescaled, etch_rad_rescaled, args['nG'], args['Lx']*1e3, loss, args['thick_max']), data, allow_pickle = True)

def main():    
    # Main: run sweep over period and losses of interest 
    for period in np.array([0.43]):
        for loss in np.array([1e-6]):
            runsweep(period, loss)

if __name__ == "__main__":
    main()