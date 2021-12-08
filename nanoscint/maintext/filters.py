import autograd.numpy as np
import autograd.scipy.signal as signal
import sys
sys.path.append('../')
import matplotlib.pyplot as plt
from scipy.signal import convolve2d as nonautograd_convolve2d

## This file contains filter functions for RCWA optimization 

def triang_to_square_tensor(Nx):
    '''
    Writes tensor to convert triangular-indexed (linear) array into 2D matrix with C4v symmetry
    This only converts a single quadrant
    '''
    idim = Nx // 2
    kdim = Nx//2*(Nx//2+1)//2
    ten = np.zeros((idim, idim, kdim))
    for i in range(idim):
        for j in range(i+1):
            for k in range(kdim):
                # if (k == j + (i*(i+1))//2 or k == i + (j*(j+1))//2):
                if (k == j + i*(i+1)//2):             
                # if (k == i + j*(j+1)//2):             
                    ten[i,j,k] = 1 
                    ten[j,i,k] = 1
    return ten

def binarization_filter(DOF, beta):
    '''
    This function assumes DOF in [0, 1]
    Steepness parameter beta 
    '''
    eps = 1e-12
    tanh_arg = np.multiply(beta, np.log(np.divide(DOF+eps, 1-DOF+eps))) 
    return np.multiply(1+np.tanh(tanh_arg), 1/2.)

def binarization_filter_hard(DOF):
    '''
    This function assumes DOF in [0, 1]
    Equivalent to beta = np.inf
    '''    
    return 0.*(DOF <= 0.5) + 1.*(DOF>0.5)

def binarization_filter_extended(DOF, beta):
    '''
    This function assumes DOF in [-\infty, +\infty]
    Steepness parameter beta 
    '''

    # Brings DOF back to [0, 1]
    DOF = np.minimum(np.maximum(DOF, 0), 1)
    tanh_arg = np.multiply(beta, np.log(np.divide(DOF, 1-DOF))) 
    return np.multiply(1+np.tanh(tanh_arg), 1/2.)    

def density_filter(DOF, rad, alpha, mode = 'constant'):
    '''
    This function assumes DOF in [0, 1]
    DOF is a 2D matrix representing full DOF (boundary conditions are used)
    Steepness parameter beta 
    rad is (circular) distance to furthest away pixel 
    alpha is steepness parameter 
    alpha = 0: single-pixel average
    alpha = +\infty : all pixels are uniformly averaged 
    '''    
    if not(rad % 2 == 0):
        raise ValueError("rad must be an EVEN integer.")
    if mode == 'constant':
        return signal.convolve(np.pad(DOF, rad//2, mode = 'constant'), gaussian_kernel(rad, alpha, False), mode = 'valid')    
    if mode == 'same':
        ''' CAREFUL '''
        ''' This method is not autograd compatible (only mode = constant is supported) '''        
        return nonautograd_convolve2d(DOF, gaussian_kernel(rad, alpha), boundary='wrap', mode='same')
    

def gaussian_kernel(size, sigma, verbose = False):
    # kernel_1D = np.linspace(-(size // 2), size // 2, size)
    kernel_1D = np.concatenate((np.arange(0, size // 2 + 1), np.arange(-(size // 2), 0)))
    # print(kernel_1D)
    kernel_1D = dnorm(kernel_1D, 0., sigma)        
    kernel_2D = np.outer(kernel_1D.T, kernel_1D.T)
    kernel_2D = np.roll(kernel_2D, kernel_2D.shape[0]//2, axis = 0)
    kernel_2D = np.roll(kernel_2D, kernel_2D.shape[1]//2, axis = 1)
    # kernel_2D *= 1.0 / kernel_2D.max()
    kernel_2D *= 1.0 / np.sum(np.sum(kernel_2D))
    
    if verbose:
        plt.imshow(kernel_2D, interpolation='none', cmap='gray')
        plt.title("Kernel ( {}X{} )".format(size, size))
        plt.colorbar()
        plt.show()    
 
    return kernel_2D

def dnorm(x, mu, sd):
    return 1 / (np.sqrt(2 * np.pi) * sd) * np.e ** (-np.power((x - mu) / sd, 2) / 2)    

def dof_to_pillar(rad, Nx, Ny, Lx, Ly, binary_flag = False):
    '''
    This function takes in a rad from DOF 
    and returns square DOF with circle  
    Nx, Ny are assumed to be even integers 
    '''        
    if not(Nx % 2 == 0 and Ny % 2 == 0):
        raise ValueError("Nx, Ny must be EVEN integers.")    
    xvec = np.linspace(-Lx/2., Lx/2., Nx)
    yvec = np.linspace(-Ly/2., Ly/2., Ny)
    [XX, YY] = np.meshgrid(xvec, yvec)    
    sigma = rad/2.355 # 2sqrt(2log(2))
    if not(binary_flag):
        DOF = np.exp(-(np.square(XX)+np.square(YY))/2./sigma**2.)
    else:
        DOF = np.zeros_like(XX)
        DOF[np.where(XX**2. + YY**2.<=rad**2.)] = 1.
    return DOF

def sin_profile(Nx, Ny, Lx, Ly):
    '''
    This function returns a perfect sin profile (given period)
    Nx, Ny are assumed to be even integers 
    '''        
    xvec = np.linspace(0, Lx, Nx)
    yvec = np.linspace(0, Ly, Ny)
    [XX, YY] = np.meshgrid(xvec, yvec)    
    DOF = np.sin(np.pi*XX/Lx)**2.*np.sin(np.pi*YY/Ly)**2.
    return DOF    

def C4v(DOFt, Nx, Ny):
    '''
    This function takes in triangular DOFt (list) 
    and returns square DOF with C4v symmetry Nx * Ny pixels
    Nx, Ny are assumed to be even integers 
    '''    
    if not(Nx % 2 == 0 and Ny % 2 == 0):
        raise ValueError("Nx, Ny must be EVEN integers.")
    # Step 1 : symmetrize matrix
    # DOF_temp = np.add(DOFt.reshape((Nx//2, Nx//2)), DOFt.reshape((Nx//2, Nx//2)).T)/2.
    # print("Debug: {} {}".format(np.shape(triang_to_square_tensor(Nx)), np.shape(DOFt)))
    DOF_temp = np.matmul(triang_to_square_tensor(Nx), DOFt)
    # DOF_temp = (DOF_temp + np.transpose(DOF_temp))/2.
   
    # Step 2 : complete C4v domain with x/y symmetries
    DOF_half = np.concatenate((np.fliplr(DOF_temp), DOF_temp), axis = 1)
    DOF = np.concatenate((np.flipud(DOF_half), DOF_half), axis = 0)
    return DOF


def inverse_C4v_sq(DOF, Nx, Ny):
    '''
    This function takes in square DOF (2d matrix) 
    and returns square DOFt 
    Nx, Ny are assumed to be even integers 
    '''    
    if not(Nx % 2 == 0 and Ny % 2 == 0):
        raise ValueError("Nx, Ny must be EVEN integers.")
    return DOF[Nx//2:, Nx//2:].reshape(((Nx//2)**2,)).flatten()

def inverse_C4v(DOF, Nx, Ny):
    '''
    This function takes in square DOF (2d matrix) 
    and returns triangular DOFt 
    Nx, Ny are assumed to be even integers 
    '''    
    if not(Nx % 2 == 0 and Ny % 2 == 0):
        raise ValueError("Nx, Ny must be EVEN integers.")
    DOFt = []
    for nx in range(int(Nx/2), Nx):
        for ny in range(int(Nx/2), nx+1):
            DOFt.append(DOF[nx, ny])
    return np.asarray(DOFt)

''' 
## FILTERS FOR 1D DOF OPTIMIZATION     
'''

def sym1d_matrix(Nx):
    mat = np.zeros((2*Nx, Nx))
    for ii in range(2*Nx):
        if ii < Nx:
            mat[ii,Nx-ii-1] = 1.
        else:
            mat[ii,ii-Nx] = 1.
    return mat

def sym1d(dof, Nx):
    # Makes 1D dof symmetric
    # return np.concatenate((np.flip(dof), dof))
    return np.matmul(sym1d_matrix(Nx), dof)

def gaussian_kernel1d(size, sigma, verbose = False):
    kernel_1D = np.concatenate((np.arange(0, size // 2 + 1), np.arange(-(size // 2), 0)))
    kernel_1D = dnorm(kernel_1D, 0., sigma)
    kernel_1D = np.roll(kernel_1D, kernel_1D.shape[0]//2, axis = 0)
    kernel_1D *= 1.0 / np.sum(kernel_1D)

    if verbose:
        plt.plot(kernel_1D)
        plt.title("Kernel ( {})".format(len(kernel_1D)))
        plt.show()    

    return kernel_1D

def density_filter1d(DOF, rad, alpha):
    '''
    This function assumes DOF in [0, 1]
    DOF is a 1D array representing full DOF (boundary conditions are used)
    rad is (circular) distance to furthest away pixel 
    alpha is steepness parameter 
    alpha = 0: single-pixel average
    alpha = +\infty : all pixels are uniformly averaged 
    '''    
    if not(rad % 2 == 0):
        raise ValueError("rad must be an EVEN integer.")
    return signal.convolve(np.pad(DOF, rad//2, mode = 'constant'), gaussian_kernel1d(rad, alpha), mode = 'valid')