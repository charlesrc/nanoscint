<div id="top"></div>

<!-- PROJECT LOGO -->
<br />

<h3 align="center">A framework for scintillation in nanophotonics</h3>

  <p align="justify">
    Code and data for paper "A framework for scintillation in nanophotonics", preprint available at https://arxiv.org/abs/2110.11492 
    If you use methods in this repository or find them useful, please cite our work as [Roques-Carmes, Charles, et al. "A general framework for scintillation in nanophotonics." arXiv preprint arXiv:2110.11492 (2021).]
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
  <li>Experimental data from electron-beam-induced scintillation from silica defects in a silicon-on-insulator photonic crystal sample </li>
  <li>Experimental data from X-ray-induced scintillation from cerium-doped yttrium aluminum garnet</li>
  <li>Rate equation model to interpret energy and current dependences in electron-beam-induced scintillation</li>    
  <li>Density functional theory (DFT) data for the spectrum and energy levels of silica defects</li>        
  <li>Numerical modeling of the two experimental configurations mentioned above, consisting in the combination of:</li>        
    <ul>
        <li> High-energy particle energy loss calculations </li> 
        <li> DFT calculations of the emitters' spectrum and energy levels </li>         
        <li> Light emission in nanostructured media </li>                 
    </ul>        
</ul>

<p align="right">(<a href="#top">back to top</a>)</p>

### Numerical methods: implementation

Our numerical methods utilize some existing libraries and packages:

* [CASINO](https://www.gel.usherbrooke.ca/casino/What.html)
* [grcwa](https://github.com/weiliangjinca/grcwa)
* [JDFTx](https://jdftx.org/)
* [NLopt](https://nlopt.readthedocs.io/en/)

<p align="right">(<a href="#top">back to top</a>)</p>


<!-- GETTING STARTED -->
## Getting Started

To get a local copy up and running follow these simple example steps.

### Prerequisites

To download required Python packages, you can use the following command
* npm
  ```sh
  npm install npm@latest -g
  ```

### Installation

1. Clone the repo
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

All inquiries should be addressed to the corresponding authors of the manuscript:

Charles Roques-Carmes - [@twitter_handle](https://twitter.com/twitter_handle) - chrc@mit.edu

Project Link: [https://github.com/github_username/repo_name](https://github.com/github_username/repo_name)

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- ACKNOWLEDGMENTS -->
## Acknowledgments

* []()
* []()
* []()

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/github_username/repo_name.svg?style=for-the-badge
[contributors-url]: https://github.com/github_username/repo_name/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/github_username/repo_name.svg?style=for-the-badge
[forks-url]: https://github.com/github_username/repo_name/network/members
[stars-shield]: https://img.shields.io/github/stars/github_username/repo_name.svg?style=for-the-badge
[stars-url]: https://github.com/github_username/repo_name/stargazers
[issues-shield]: https://img.shields.io/github/issues/github_username/repo_name.svg?style=for-the-badge
[issues-url]: https://github.com/github_username/repo_name/issues
[license-shield]: https://img.shields.io/github/license/github_username/repo_name.svg?style=for-the-badge
[license-url]: https://github.com/github_username/repo_name/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/linkedin_username
[product-screenshot]: images/screenshot.png
