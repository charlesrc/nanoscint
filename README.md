<div id="top"></div>

<!-- PROJECT LOGO -->
<br />

<h3 align="center">A framework for scintillation in nanophotonics</h3>

  <p align="justify">
    Code and data for paper "A framework for scintillation in nanophotonics", preprint available on <a href="https://arxiv.org/abs/2110.11492">arXiv</a>. If you use methods in this repository or find them useful, please cite our work as [Roques-Carmes, Charles, et al. "A general framework for scintillation in nanophotonics." arXiv preprint arXiv:2110.11492 (2021).]
    <br />
  </p>
</div>

<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About</a>
      <ul>
        <li><a href="#built-with">Numerical methods</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>


<!-- ABOUT THE PROJECT -->
## About

This repository contains codes and data published in the above-mentioned manuscript. The methods presented here can be used to model, analyze, and optimize nanophotonic scintillators. Specifically, this repository contains: 

<ul>
  <li>Numerical modeling of the two experimental configurations mentioned above, consisting in the combination of:</li>        
    <ul>
        <li> High-energy particle energy loss calculations </li> 
        <li> DFT calculations of the emitters' spectrum and energy levels </li>         
        <li> Light emission in nanostructured media </li>                 
    </ul>        
  <li>Experimental data from electron-beam-induced scintillation from silica defects in a silicon-on-insulator photonic crystal sample </li>
  <li>Experimental data from X-ray-induced scintillation from cerium-doped yttrium aluminum garnet</li>
  <li>Rate equation model to interpret energy and current dependences in electron-beam-induced scintillation</li>    
  <li>Density functional theory (DFT) data for the spectrum and energy levels of silica defects</li>        
  <li>Methods for topology optimization of nanophotonic scintillators with rigorous coupled-wave analysis (RCWA)</li>          
</ul>

<p align="right">(<a href="#top">back to top</a>)</p>

### Numerical methods: implementation

Our numerical methods utilize some existing libraries and packages:

* [grcwa](https://github.com/weiliangjinca/grcwa) for photonics simulation and automatic differentiation with RCWA.
* [CASINO](https://www.gel.usherbrooke.ca/casino/What.html) for electron-beam energy loss calculations. Examples of output data from CASINO used in our work can be found [here](nanoscint/maintext/Figure2/casinores/). 
* [JDFTx](https://jdftx.org/) for DFT calculations. Examples of output data from JDTFx used in our work can be found [here](nanoscint/maintext/DFT Setup and Analysis/). 
* [NLopt](https://nlopt.readthedocs.io/en/latest/) for optimization.

If you use some of those methods, please cite them as well, as appropriate. 

<p align="right">(<a href="#top">back to top</a>)</p>


<!-- GETTING STARTED -->
## Getting Started

To get a local copy up and running follow these simple example steps.

### Python Prerequisites

Code was written in Python 3.9.6. To download required Python packages, you can use the following command (preferably in a virtual environment):
* Python Requirements
  ```sh
  pip install -r requirements.txt
  ```

### Julia Prerequisites

Code was written in Julia 1.6.1. To download required Julia packages, you can use the following command (preferably in a virtual environment):
* Julia Requirements
  ```sh
  Pkg.add(["PyCall", "PyPlot", "Peaks", "Statistics", "LinearAlgebra", "SparseArrays", "GSL", "DifferentialEquations", "Sundials", "PyCall", "LsqFit", "Printf", "Interpolations", "DelimitedFiles", "JLD2"])
  ```

### Installation

Clone the repository
   ```sh
   git clone https://github.com/charlesrc/nanoscint.git
   ```

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#top">back to top</a>)</p>


<!-- CONTACT -->
## Contact

All questions, inquiries, or suggestions should be addressed to the corresponding authors of the manuscript:

* Charles Roques-Carmes - [@personal_website](https://chrc.scripts.mit.edu) - chrc@mit.edu
* Nicholas Rivera - [@personal_website](http://nrivera.scripts.mit.edu/nhr/) - nrivera@mit.edu 
* Project Link: [https://github.com/charlesrc/nanoscint](https://github.com/charlesrc/nanoscint)

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- ACKNOWLEDGMENTS -->
## Acknowledgments

* This material is based upon work supported in part by the U.S. Army Research Laboratory and the U.S. Army Research Office through the Institute for Soldier Nanotechnologies, under contract number W911NF-18–2–0048. 
* This material is also in part based upon work supported by the Air Force Office of Scientific Research under the award number FA9550-20-1-0115, as well as in part supported by the Air Force Office of Scientific Research under the award number FA9550-21-1-0299.
* This work was performed in part on the Raith VELION FIB-SEM in the MIT.nano Characterization Facilities (Award: DMR-2117609)

<p align="right">(<a href="#top">back to top</a>)</p>
