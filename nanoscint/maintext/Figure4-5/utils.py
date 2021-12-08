import nlopt
import autograd.numpy as np
import sys
import os

def rescaled(dof, dof_min, dof_max):
    return dof_min + dof*(dof_max-dof_min)