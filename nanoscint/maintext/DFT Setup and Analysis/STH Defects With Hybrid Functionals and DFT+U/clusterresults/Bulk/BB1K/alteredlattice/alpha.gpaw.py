#! /usr/bin/env python3

import numpy as np
from ase import Atoms
import ase.units as units
from gpaw import GPAW, PW
from gpaw.hyperfine import hyperfine_parameters
from ase.io import write
from math import sqrt
from gpaw import GPAW, PW, FermiDirac

from ase.units import Bohr

alpha = Atoms('OOOOOOSiSiSi', scaled_positions=[(.43, .266, .214), (.732, .1462, .54733), (0.8538, 0.587, 0.880667), (0.2668, 0.413, 0.786), 
(0.1462, 0.7332, 0.452667), (0.587, 0.8538, 0.11933), (0.4699, 0, 0.3333), (0, 0.4688, 0.6666667), (0.5301, 0.5301, 0)],
 cell = [Bohr*np.array([4.914, -8.51129, 0]), Bohr*np.array([4.914, 8.5113, 0]), Bohr*np.array([0, 0, 10.82])], pbc=True)

alpha_p = Atoms('OOOOOOSiSiSi', scaled_positions=[(.563, .2668, .2314), (.732, .1462, .54733), (0.8538, 0.587, 0.880667), (0.2668, 0.413, 0.786),
(0.1462, 0.7332, 0.452667), (0.587, 0.8538, 0.11933), (0.4699, 0, 0.3333), (0, 0.4688, 0.6666667), (0.5301, 0.5301, 0)],
 cell = [Bohr*np.array([4.914, -8.51129, 0]), Bohr*np.array([4.914, 8.5113, 0]), Bohr*np.array([0, 0, 10.82])], pbc=True, magmoms=[1, 0, 0, 0, 0, 0, 0, 0, 0])


bn22bc = Atoms('CNBBBNNN', scaled_positions=[(0, 0, 0), (1/3, -1/3, 0), (0.5, 0, 0), (0, 0.5, 0), (0.5, 0.5, 0), (1/3+0.5, -1/3, 0), (1/3,-1/3+0.5, 0), (1/3+0.5, -1/3+0.5, 0)], cell=[[1.42*sqrt(3)*2, 0, 0], [-1.42*sqrt(3), 1.42*3, 0], [0, 0, 10]], pbc = [1, 1, 0])

write('bn22bc.png', bn22bc)

alpha_p.calc = GPAW(mode=PW(500), xc='PBE', kpts=(1, 1, 1), occupations=FermiDirac(0.01, fixmagmom=True), txt=None, charge=+1)

#e = alpha_p.get_potential_energy()
print(alpha_p.get_magnetic_moments())
print(alpha_p.calc.get_spin_polarized())


write('ASE_Alpha.png', alpha)
write('ASE_Alphap.png', alpha_p)

#h = Atoms('H', magmoms=[1])
#h.center(vacuum=3)
#h.calc = GPAW(mode=PW(1000), txt=None)
#e = h.get_potential_energy()
#A = hyperfine_parameters(h.calc)[0] * 5.586
#a = np.trace(A) / 3
#frequency = a * units._e / units._hplanck  # Hz
#wavelength = units._c / frequency  # meters
#print(f'{frequency*1/1000000} MHz')

#print(f'{wavelength * 100:.1f} cm')

#lattice 
#4.914  4.914 0.00000000000000 \
#-8.51129766839346  8.51129766839346  0.00000000000000 \
#0.00000000000000   0.00000000000000   10.81199999997936

#ion O 0.56300000000000 0.26680000000000  0.231400000000000  0




#ion O 0.43300000000000 0.26680000000000  0.21400000000000  0
#ion O 0.73320000000000 0.1462000000000  0.54733300000000 0
#ion O 0.85380000000000   0.5870000000000   0.88066700000000 0
#ion O 0.2668000000000   0.4130000000000  0.7860000000000 0
#ion O 0.1462000000000  0.73320000000000 0.45266700000000 0
#ion O 0.5870000000000 0.85380000000000  0.11933300000000 0
#ion Si 0.46990000000000 0.0000000000000 0.33333300000000 0
#ion Si  0.00000000000000 0.46990000000000 0.66666700000000 0
#ion Si  0.53010000000000  0.5301000000000 0.00000000000000  0

