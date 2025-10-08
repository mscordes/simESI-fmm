![Cover Image](/assets/cover_image.jpg)
# simESI-fmm
simESI-fmm (simulations of ESI with FMM) is a high performance version of the original "simESI" (https://github.com/mscordes/simESI). As with the original simESI release, simESI-fmm is a program utilizing GROMACS for simulating electrospray ionization (ESI) of proteins in ammonium acetate containing droplets to form protonated or deprotonated protein ions. simESI-fmm enables this by allowing proton transfer reactions between discrete amino acids, water, Grotthuss diffuse H₃O⁺ and OH⁻, ammonium (NH₄⁺), ammonia (NH₃), acetate (CH₃COO⁻), and acetic acid (CH₃COOH). Additional models are included to enable modelling of ambient conditions. simESI-fmm handles the simulation from end-to-end including preprocessing of inputted protein coordinate file, droplet formation & equilibration, and running of the simulation. The simulation can be readily modified in many ways from command line (see below). 

## Changes from the original simESI 
simESI-fmm was designed from the ground up to enable simulation at scale, i.e., larger droplets, larger proteins, and much better handling of protein complexes. Relative to simESI, the numerical calculations surrounding proton transfer are much quicker. Internal simESI calculations occupy a negligible portion of the runtime relative to the gmx calls, particularly at scale.

Additionally, simESI-fmm optionally supports the use of the fast multipole method (FMM) for calculation of non-bonded forces as the name suggests which can significantly speed the MD simulations for larger systems. 

## Citations
* Cordes, M. S. & Gallagher, E. S. Molecular Dynamics Simulations of Native Protein Charging in Electrosprayed Droplets with Experimentally Relevant Compositions. *J. Am. Chem. Soc.* 147, 15066-15076 (2025). https://doi.org/10.1021/jacs.4c17382
* Kohnke, B., Kutzner, C. & Grubmüller, H. A GPU-Accelerated Fast Multipole Method for GROMACS: Performance and Accuracy. *J. Chem. Theory Comput.* 16, 6938-6949 (2020). https://doi.org/10.1021/acs.jctc.0c00744

## Dependencies
* GROMACS, tested on GROMACS 2022.3, but should work on any fairly recent version
 ** Optionally, if using GROMACS with FMM see https://grubmueller.pages.mpcdf.de/docs-gromacs-fmm-constantph/docs/install_guide/
* Highly recommended that the python module PROPKA3 is callable from command line for residue specific pKa values, else default pKa's will be used.
 ** Should be able to call ```python -m propka``` from the terminal where simESI-fmm is being used.

## Installation 
After downloading, navigate to the ```simESI-fmm``` parent directory and run the following commands as simESI-fmm must be compiled prior to running. simESI-fmm has been extensively tested on linux and windows but should work on mac as well.

For windows,
```mkdir build``` <br />
```cd build``` <br />
```cmake ..``` <br />
```cmake --build . --config Release``` <br />

For linux,
```mkdir build``` <br />
```cd build``` <br />
```cmake -DCMAKE_BUILD_TYPE=Release ..``` <br />
```cmake --build .``` <br />

If compilation is successfull, this should create a ```simESI``` executable in the ```simESI-fmm``` parent directory.

## Running 
To run simESI-fmm with the default example, simply navigate to the ```simESI-fmm``` parent directory and call ```./simESI``` This will create a trial output directory like ```ubq_1``` in ```outputFiles/``` with the default ```ubq.pdb``` structure included in ```coordinateFiles/```. Each subsequent call will count up to ```ubq_2```, ```ubq_3```, and so on.

Within each trial output directory there are two subdirectories, ```data``` and ```simulation```. ```data``` contains ```.txt``` files with all the high level information to expedite data analysis including information like protein charge, system charge, number of each molecule, and execution times. ```simulation``` contains all the files from the simulation itself, namely coordinate (```.gro```) files from each step. 

To adjust simulation parameters, simESI uses a number of command line args detailed below to modify the simulation and how the system/droplet are defined.  

## Run Parameters
To see a description of run parameters call ```./simESI -h```, or open ```inputFiles/default.txt```. simESI-fmm uses command line arguments and/or an input file system to determine run parameters.

### Command Line Arguments
As an example of how command line arguments can be used, if you would like to run with a different protein, place the associated ```.pdb``` (must be ```.pdb```!) of the protein files in ```coordinateFiles/``` and call ```./simESI -pdb``` with the new ```.pdb``` name. 

### Input File System
simESI-fmm can also extract run parameters from input files stored in the ```inputFiles/``` subdirectory. The default parameters are stored within ```inputFiles/default.txt```, which is also an example for how input files should be formatted (NOTE: always keep a ```default.txt``` in ```inputFiles/```). To modify defaults, simply modify this file.

To make a new input file, make a copy of ```default.txt```, name it something else, and alter the parameters as you see fit. To call the new input file, run ```./simESI -f``` with the name of the new input file. Input file format is ```arg```, then ```=``` sign, then ```value```. Anything after a # character is a comment and will be disregarded. Note that by default, command line arguments will override input file parameters.

## Data Analysis
simESI-fmm continually saves data associated from each step in the run in the data subdirectory. Included in the ```plotting/``` subdirectory is ```composition.py``` that demonstrates how data from runs like droplet composition can be analyzed. Additionally, ```movie.py``` is included that produces frames from individual simESI-fmm ```.gro``` files which can be stitched together to create a movie.

* ```composition.py``` This script plots the evolution of droplet composition in addition to protein charge, and ratio of net droplet charge to the droplets Rayleigh limit. To choose the trial to plot, input the trial directory in ```outputFiles/``` to ```python composition.py -d ubq_1 -m 8.6``` (ie. ```ubq_1``` for a default run with ubiquitin). ```--m``` in this case is the protein mass (in kDa) which is required to compute the Rayleigh limit.

<p align="center">Example output from composition.py using ubiquitin as an example.</p>

<p align="center">
  <img src="/assets/sysinfo_ex_ubq.png" width="500">
</p>

## Notes
* **simESI-fmm cannot use conventional trajectories between steps!** This is an unfortunate byproduct of continually changing the number and types of molecules due to proton transfers. For posterity, simESI-fmm saves the coordinate file outputted after each (4 ps) step. A script is included in this package, ```plotting/movie.py``` that create frames from generated ```.gro``` files in pymol which can be stitched together in order to make a movie.
* **simESI-fmm is memory intensive!** Another unfortunate byproduct of the above point is that given the number of steps over an entire simulated run, many coordinate files are saved. A complete run with ubiquitin which is a very small protein will be a couple of GB's, larger proteins will be even more.
* *Ammonium acetate concentration has a significant effect on the final ion charge states produced by simESI-fmm. As a rule of thumb, higher ammonium acetate concentration produces lower charge states and vice versa. If charge states are different consider modifying the default initial, droplet ammonium acetate concentration of 0.25 M. For example, to increase to 0.50 M you would call ```./simESI -amac 0.50```.*
