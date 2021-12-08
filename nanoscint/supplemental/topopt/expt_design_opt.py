# Written for Python 3.9.6
# Packages 
import nlopt # For non-linear optimization. See download instructions at https://nlopt.readthedocs.io/en/latest/
import autograd.numpy as np # Important! See download instructions at https://github.com/HIPS/autograd
from autograd import grad 
import sys
import matplotlib.pyplot as plt 

sys.path.append('../../maintext/') 
import filters # Topology optimization filter functions 

import grcwa # Download instructions at https://github.com/weiliangjinca/grcwa
grcwa.set_backend('autograd')  # Important to make grcwa work with autograd 

from scipy.io import savemat

args = {}
params = {}

## Define variables 
args['nG'] = 51 # Always check convergence wrt nG

# Lattice vector (unit length is 1 um)
args['Lx'] = 0.430
args['Ly'] = args['Lx']
args['lam0'] = 0.670
args['freq0'] = 1./args['lam0']
args['L1'] = [args['Lx'],0.]
args['L2'] = [0.,args['Ly']]

# Wavelength vector 
wl_vec = np.linspace(0.4,0.8,400)
        
# discretization for patterned layers
# Has to be even integers with current C4v implementation
args['Nx'] = 50
args['Ny'] = args['Nx']

## Multilayer setting. Bottom/Top are half infinite spaces, middle is for unpatterned middle layers. 
eps_SiO2 = 2.13353 +1e-5*1j# IP-Dip
eps_Si = 14.63 + 0.0823*1j
eps_vacuum = 1.
eps_bot = 1.
eps_top = 1.

# Set for initialization 
# For now assumes typical multilayer setup 
# Define min and max parameters
args['thick_max'] = 0.5
args['etch_depth'] = 0.45
thick_SiO2 = 0.3

theta = 0.
phi = 0.

# Filter parameters
args['rad'] = 8 # max. radius of density filter 
args['alpha'] = 5.0 # std dev of density filter 
args['beta'] = 1000.0 # scaling parameter of binarization filter 

eps_Si_data = np.loadtxt(open('Si (Silicon) - Palik.txt', "rb"), delimiter=",", skiprows=0)

def eps_Si_func(wl):
    ''' Permittivity of silicon function '''
    ind = np.where(np.abs(eps_Si_data[:,0]-wl) == np.min(np.abs(eps_Si_data[:,0]-wl)))       
    return (eps_Si_data[ind,1] + 1j * eps_Si_data[ind,2])[0][0]

def eps_SiO2_func(wl):
    ''' Permittivity of silica function '''
    width = 0.078
    f = 1/wl
    f0 = 1/0.670
    # Fake red peak small losses 
    epsi = 2.45e-5 * np.exp(-np.divide(np.square(f-f0),2.*width**2.))
    return 2.13353 + 1j * epsi    

def rescaled(dof, dof_min, dof_max):
    return dof_min + dof*(dof_max-dof_min)

def Veff(freqangs, dof, ref_flag = False):
    freq = freqangs[0]
    theta = freqangs[1]
    phi = freqangs[2]

    eps_Si = eps_Si_func(1/freq)

    ## All DOFs are assumed to be in [0, 1]
    obj = grcwa.obj(args['nG'], args['L1'], args['L2'], freq, theta, phi, verbose = 0)      
    obj.Add_LayerUniform(0., eps_vacuum)
    obj.Add_LayerGrid(args['etch_depth'], args['Nx'], args['Ny'])
    obj.Add_LayerUniform(args['thick_max']-args['etch_depth'], eps_Si)
    obj.Add_LayerUniform(thick_SiO2, eps_SiO2)
    obj.Add_LayerUniform(100., eps_Si)
    obj.Add_LayerUniform(0., eps_vacuum)

    obj.Init_Setup(Gmethod = 0)
    
    # Updates DOF 
    flattened_dof = filters.C4v(dof, args['Nx'], args['Ny'])
    if not(ref_flag):
        flattened_dof = filters.density_filter(flattened_dof, args['rad'], args['alpha'])
        flattened_dof = filters.binarization_filter(flattened_dof, args['beta'])    
    flattened_dof = np.array(flattened_dof).flatten()*(eps_Si - eps_vacuum) + eps_vacuum # Turns 0/1 DOF into epsilon values        
    obj.GridLayer_geteps(flattened_dof)

    planewave = {'p_amp': 0, 's_amp': 1, 'p_phase': 0, 's_phase': 0}
    obj.MakeExcitationPlanewave(planewave['p_amp'], planewave['p_phase'], planewave['s_amp'], planewave['s_phase'], order = 0)
    
    # obj.Print()
    SiO2_layer = 3
    dN = 1/args['Nx']/args['Ny']
    M0 = grcwa.fft_funs.get_conv(dN, np.ones((args['Nx'], args['Ny'])), obj.G)
    vol1 = obj.Volume_integral(SiO2_layer, M0, M0, M0, normalize=1)

    res = np.abs(vol1) 

    return res

gradVeff = grad(Veff,1)
freqangs0 = [args['freq0'], theta, phi]
iter_count = 0
def cost(dof, g):
    global iter_count
    global best_cost
    iter_count+=1
    if g.size > 0:
        g[:] = gradVeff(freqangs0,dof)
        print(g)
    res = Veff(freqangs0, dof)
    if iter_count % 10 == 0 :
        data = {'params': params, 'args': args, 'dof': dof, 'freq_vec': args['freq0'], 'theta_vec': theta, 'phi_vec': phi, 'cost':res}
        np.save('res/' + folder_name + '/dof_{}'.format(iter_count), data, allow_pickle = True)        
    if res > best_cost:
        best_cost = res
        data = {'params': params, 'args': args, 'dof': dof, 'freq_vec': args['freq0'], 'theta_vec': theta, 'phi_vec': phi, 'cost':res}
        np.save('res/' + folder_name + '/dof_{}'.format("best"), data, allow_pickle = True)        
    print("Obj. fun. iter {}, cost = {}, enhanc. = {}".format(iter_count, res, res/ref))
    return res

folder_name = 'dof_opt_4'

iter_count = 0             

nDOF = args['Nx']//2 * (args['Nx']//2 +1)//2 # just thicknesses 

opt = nlopt.opt(nlopt.LD_MMA, nDOF)
print("Number of DOFs: {}".format(nDOF))

tol = 1e-8

file_init = "" # No initial data (otherwise change to name of dof file )
print_flag = False # Calculates and prints performance of file_init data or not 

if file_init != "":
    file_reader = np.load("res/" + folder_name + "/" + file_init, allow_pickle = True)
    dof_init = file_reader.item().get('dof')
    best_cost = file_reader.item().get('cost')
else:
    dof_init = np.random.rand(nDOF,)
    best_cost = 0.

opt.set_max_objective(cost)

# Varying thicknesses, varying DOFs 
opt.set_upper_bounds(np.ones(nDOF,))
opt.set_lower_bounds(np.zeros(nDOF,))

ref = Veff(freqangs0, np.ones((nDOF,)), ref_flag = True)
print("REFERENCE ABSORPTION FOR UNPATTERNED = {}".format(ref))

if not(print_flag): 
    # Runs topology optimization 
    opt.set_xtol_rel(tol)
    dof_sol = opt.optimize(dof_init)
    optf = opt.last_optimum_value()
    print("optimum at ", dof_sol)
    print("optimal value = ", optf)
    print("result code = ", opt.last_optimize_result())
    data = {'params': params, 'args': args, 'dof': dof_sol, 'freq_vec': args['freq0'], 'theta_vec': theta, 'phi_vec': phi, 'cost':optf}
    np.save('res/' + folder_name + '/dof_{}_best'.format(iter_count), data, allow_pickle = True)
else:
    # Binarizes dof and calculates performance of design 
    dof_sol = dof_init
    print("Optimal dof {}".format(dof_sol))
    dof_sol_bin = filters.C4v(dof_sol, args['Nx'], args['Ny'])
    dof_sol_bin = filters.density_filter(dof_sol_bin, args['rad'], args['alpha'])
    dof_sol_bin = filters.binarization_filter(dof_sol_bin, args['beta'])    
    dof_sol_bin = filters.binarization_filter_hard(dof_sol_bin)    
    dof_sol_bin = filters.inverse_C4v(dof_sol_bin, args['Nx'], args['Ny'])
    print("Binarized dof {}".format(dof_sol_bin))

print("FINAL ABSORPTION FOR PATTERNED = {}".format(Veff(freqangs0, dof_sol)))
print("FINAL ABSORPTION FOR PATTERNED, BINARIZED = {}".format(Veff(freqangs0, dof_sol_bin, ref_flag = True)))

# Plot DOFs and binarized DOFs
plt.imshow(filters.C4v(dof_sol, args['Nx'], args['Ny']))
plt.savefig('dof_si.png')

plt.imshow(filters.C4v(dof_sol_bin, args['Nx'], args['Ny']))
plt.savefig('dof_si_bin.png')

matlab_dict = {"Nx": float(args["Nx"]), "Ny": float(args["Ny"]), "Lx": args['Lx'], "Ly": args['Ly'], "thick_max": args["thick_max"], "etch_depth": args["etch_depth"],"dof": filters.C4v(dof_sol_bin, args['Nx'], args['Ny'])}
savemat('res/' + folder_name + '/dof_si_bin.mat', matlab_dict)

Veff_vec= []
eps_SiO2_vec = []
for wl in wl_vec:
    freqangs0 = [1/wl, theta, phi]
    res = Veff(freqangs0, dof_sol_bin)/wl
    eps_SiO2_vec.append(eps_SiO2_func(wl))    
    print("Computing wl = {} nm".format(wl*1e3))
    Veff_vec.append(res)

# Plot absorption spectrum
plt.figure()
plt.plot(wl_vec, Veff_vec)