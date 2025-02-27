![Cover Image](/assets/cover_image.jpg)
# simESI
NOTE: This program is under active development.
simESI-fmm (simulations of ESI with FMM) is a high performance version of the original "simESI" (https://github.com/mscordes/simESI). As with the original simESI release, simESI-fmm is a program utilizing GROMACS for simulating electrospray ionization (ESI) of proteins in ammonium acetate containing droplets to form protonated or deprotonated protein ions. simESI-fmm enables this by allowing proton transfer reactions between discrete amino acids, water, Grotthuss diffuse H₃O⁺ and OH⁻, ammonium (NH₄⁺), ammonia (NH₃), acetate (CH₃COO⁻), and acetic acid (CH₃COOH). Additional models are included to enable modelling of ambient conditions. simESI-fmm handles the simulation from end-to-end including preprocessing of inputted protein coordinate file, droplet formation & equilibriation, and running of the simulation. The simulation can be readily modified in many ways from command line (see below). 

## Changes from the original simESI distribution
simESI-fmm was designed from the ground up to enable simulation at scale, i.e., larger droplets, larger proteins, and much better handling of protein complexes. Relative to simESI, the numerical calculations surrounding proton transfer are much quicker. This is due to being entirely written in C++, more efficient clustering, and changes to droplet formation that can reduce the number of simulated waters. With most systems, proton transfer calculations will typically occupy a negligible contribution to total compute time (<5%). 

Additionally, simESI-fmm optionally supports the use of the fast multipole method (FMM) for calculation of non-bonded forces as the name suggests (if using FMM, please use the FMM citation below), which can significantly speed up the MD portion of the simulation given that droplets have net charge (and therefore cannot use PME). 

## Citations
If using simESI, please cite the following paper(s).

Cordes, M.S.; Gallagher, E.S. Molecular Dynamics Simulations of Native Protein Charging in Electrosprayed Droplets with Experimentally Relevant Compositions. 2025. Manuscript under review. doi: https://doi.org/10.26434/chemrxiv-2024-smnj3

Kohnke, B. Kutzner; C. & Grubmüller, H. J. Chem. Theory Comput. 16, 6938-6949 (2020)

## Installation and Running
simESI-fmm has been extensively tested on linux and to a lesser extent windows. simESI-fmm only uses the C++20 standard library, no external dependencies! simESI-fmm is untested on mac, but may work given linux support.

### Dependencies
* GROMACS, tested on GROMACS 2022.3, but should work on any fairly recent version
 ** Optionally, if using GROMACS with FMM see https://grubmueller.pages.mpcdf.de/docs-gromacs-fmm-constantph/docs/install_guide/
* For linux/mac, GNU (gcc/libstdc++) 11.3+
* For windows, MSVC(2022)
* cmake (version 3.20+)

### Installation
After downloading, navigate to the ```simESI-fmm``` parent directory and run the following commands as simESI-fmm must be compiled prior to running.

```mkdir bin``` <br />
```cd bin``` <br />
```cmake -DCMAKE_BUILD_TYPE=Release ..``` <br />
```cmake --build . --config Release``` <br />

If compilation is successfull, this should create a ```simESI-fmm``` executable in the ```simESI-fmm``` parent directory.

### Running simESI
simESI-fmm is a command line program. To run simESI-fmm with the default example, simply navigate to the ```simESI-fmm``` parent directory and call ```./simESI-fmm -pdb ubq.pdb```. This will start a droplet simulation with the default protein structure, ubiquitin (PDB: 1UBQ). simESI-fmm is designed to take an input ```.pdb``` file of the protein (must be ```.pdb```!) from the ```simESI-fmm/coordinateFiles``` subdirectory. 

This will create an output directory ```simESI-fmm/outputFiles``` if ```simESI-fmm/outputFiles``` does not exist. simESI-fmm then creates a numbered subdirectory using the name of the inputted ```.pdb``` in ```simESI-fmm/outputFiles```, for example ```simESI-fmm/outputFiles/ubq_1``` if using the default ```ubq.pdb``` structure. If you were to run ```./simESI-fmm -pdb ubq.pdb``` again, the new output directory would be ```simESI-fmm/outputFiles/ubq_2``` with each additional run counting up. Within each output dir (i.e. ```simESI-fmm/outputFiles/ubq_1```) there are two subdirectories, ```data``` and ```simulation```. ```data``` contains ```.txt``` files with all the high level information to expedite data analysis including information like protein charge, system charge, number of each molecule, and execution times at each timestep. ```simulation``` contains all the files from the simulation itself, namely coordinate files from each timestep, as well as files assocaiated with droplet creation/equilibriation. 

To adjust simulation parameters, simESI uses a number of command line args detailed below to modify the simulation and how the system/droplet are defined.  

## Notes
* **simESI-fmm cannot use conventional trajectories between steps!** This is an unfortunate byproduct of continually changing the number and types of molecules due to proton transfers. For posterity, simESI-fmm saves the coordinate (```.gro```) file outputted after each (4 ps) step, and trajectory (```.xtc```) file after every 10 steps to allow users to see how the simulation run is evolving. A script is included in this package, ```simESI-fmm/dataAnalysis/movie.py``` that create frames from generated ```.gro``` files in pymol which can be stitched together in order to make a movie.
* **simESI-fmm is memory intensive!** Another unfortunate byproduct of the above point is that given the number of steps over an entire simulated run, many coordinate files are saved. A complete run with ubiquitin which is a very small protein will be a couple of GB's, larger proteins will be even more.
* In terms of forcefield choice, **simESI-fmm is only capable of running simulations with the ```charmm36``` forcefield**. Additionally, given the number of non-standard molecules, simESI-fmm must be run with the  ```charmm36``` version contained with ```simESI-fmm/inputFiles```.
* *Ammonium acetate concentration has a significant effect on the final ion charge states produced by simESI-fmm. If charge state distributions produced from simESI-fmm are different then what would be expected experimentally, consider modifying the initial, droplet ammonium acetate concentration via the command line args detailed below. As a rule of thumb, higher ammonium acetate concentration produces lower charge states and vice versa.*

## Command Line Arguments
* ```-pdb ``` *(Required)* Input ```.pdb``` filename, must be present in the ```simESI-fmm/coordinateFiles``` directory. Must be supplied. As an example, can set to ```ubq.pdb``` that is included in simESI-fmm by default.
  * **Warning** *simESI-fmm requires that a valid ```.pka``` from ```propka``` can be generated. You should be able to run ```propka``` without error, for the inputted ```.pdb```. Very often, ```propka``` will miss the C-termini which will break the simulation. To fix this, the final oxygen, for every C-terminus in the inputted ```.pdb```, must have the atom name ```OXT```.*
  * **Warning** *The inputted ```.pdb``` must be valid! Many stock ```.pdb``` are missing atoms or improperly defined so you'll need to fix them prior to using simESI-fmm! If you can do ```gmx mdrun``` and ```python -m propka``` with a given structure outside of simESI-fmm and recieve no errors, you should be good to go.*
  * **Warning** *All HETATMS including solvents or ligands will be automatically deleted.*
  * **Note** *The inputted ```.pdb``` can be a protein homo- or hetero-complex.*
  
### Droplet Definition Arguments
* ```-droplet_size``` *(Optional, float)* **Note** ```simESI-fmm``` handles droplet formation differently then original ```simESI``` choosing instead to seed an 'envelope' of water around the protein, rather than a large enough sphere to fit around the protein. This is particulary advantageous when simulating non-globular proteins that would require very large spheres for solvation. Default envelope size in nm is 1.10. 
* ```-esi_mode``` *(Optional)* Choose to simulate a positive-mode ESI with ```pos``` or negative-mode ESI with ```neg```. Defaults to ```pos```.
* ```-amace_conc``` *(Optional, float)* Initial molar concentration of ammonium acetate to seed droplet with. Default is 0.25.

### System Definition Arguments
* ```-time``` *(Optional, float)* Time in ns in which to end the simulation. Default of 12.5 ns. This works well for ubiquitin, but for larger proteins will need to increase. 
* ```-water_cutoff``` *(Optional, int)* Once cutoff reached, ramp droplet temperature to evaporate the last handful waters. Will calculate based on protein mass default. For ubiquitin this yields 30 waters. 
* ```-atm``` *(Optional)* Choose to simulate droplet in fully atomistic atmosphere of 79% N₂ and 21% O₂ at 300K if ```yes```, or hard vacuum if ```no```. Default of ```yes```.
* ```-init_temp``` *(Optional, float)* Initial droplet temperature in K. Higher temperature = faster evaporation and much faster overall simulation, but risk denaturing protein. Default of 370.
* ```-final_temp``` *(Optional, float)* Final droplet temperature in K once ```-water_cutoff``` reached. Default of 450. For no temperature ramp set to 370.
* ```-box_size``` *(Optional, float)* Size of the simulation box (ie., how much atmosphere NOT droplet size). Default in nm of 75. This will form a fairly 75 x 75 x 75 nm box

### Compute Arguments
* ```-gpu``` *(Optional)* Choose to use GPU acceleration if ```yes``` which is *highly* recommended to the length of the simulations. Can run purely on CPU if ```no```. 
 * **Warning** *simESI-fmm assumes GPU accessibility with the default of ```yes```.*
* ```-hpc``` *(Optional)* Choose to run with more ```gmx mdrun``` params. This assumes 1 node and GPU acceleration. If using different resources/setup, you may need to manually modify some ```Popen``` calls within simESI-fmm related to ```gmx mdrun``` to reflect the resources you are pulling. Default of ```no```.
* ```-gmx_env``` *(Optional)* Path to gmx executable. If calling from ```gmx``` from command line, leave as default. If you want to use a specific build of GROMACS (say an FMM enabled build), set arg to the path of the gmx executable (i.e., ```"/users/foo/gromacs-fmm/build/bin/gmx"```). Default of ```gmx```.
* ```-fmm``` *(Optional)* Choose to use FMM enabled simESI-fmm assuming your ```gmx``` build ```-gmx_env``` is set to an FMM enabled GROMACS build. Default of ```no```.
 * **Note** FMM is much slower for smaller proteins, only recommend for large (maybe 50-100 kDa+?). With very large proteins though, FMM is *much* quicker. 

### Run Continuation Arguments
*If runs fail for whatever reason, and you need to restart from specific step in a specific trial run, use the args below to continue the run. If starting from the failed step does not work, consider starting well before it.*
 * **Note** Also, you must supply the ```.pdb ``` file to  ```-pdb ``` that was used when the dir was created in addition to ```-dir``` and ```-step```!* 
* ```-dir_cont``` *(Optional)* Name of subdirectory in ```simESI-fmm/outputFiles``` to continue from. For example, if using the default ```ubq``` this could be ```ubq_1```.
* ```-step_cont``` *(Optional)* Step number to continue from (ie. the step when the run failed). simESI-fmm will generate a fresh ```.top``` file, so you can restart from any given steps ```.gro``` file. **Warning** *If restarting from a step, data will automatically append to the end of the data files in the ```data``` subdirectory (see below), which may compromise data continuity. It's on the user to modify those files accordingly.*

## Data Analysis
simESI-fmm continually saves data associated from each step in the run in the data subdirectory (ie. for a default ubiquitin run this would be in, ```simESI-fmm/outputFiles/ubq_1/data```, see the **Running simESI-fmm** section for more details). Included in ```simESI-fmm/dataAnalysis``` subdirectory are two example scripts that demonstrate how the high level data stored in the ```data``` directory can be analyzed detailed below. Additionally, a third file ```movie.py``` is included that produces frames from individual simESI-fmm ```.gro``` files which can be stitched together to create a movie.

1. ```system_info.py``` This script plots the evolution of droplet composition in addition to protein charge, and ratio of net droplet charge to the droplets Rayleigh limit. To choose the trial to plot, input the trial directory in ```simESI-fmm/outputFiles``` to ```--dir``` (ie. ```ubq_1``` for a default run with ubiquitin). Additionally, computing the Rayleigh limit requires protein mass which is given (in kDa) to the ```--mass``` argument (ie. 8.6 for ubiquitin).

<p align="center">Example output from system_info.py using ubiquitin as an example.</p>
s
<p align="center">
  <img src="/assets/sysinfo_ex_ubq.png" width="500">
</p>

2. ```runtime.py``` This script is useful for analyzing the performance of simESI-fmm during a given trial. To choose the trial to plot, input the trials directory in ```simESI-fmm/outputFiles``` to ```--dir``` (ie. ```ubq_1``` for a default run with ubiquitin). This creates four stacked plots including the execution times from each timestep of the simESI-fmm exchange script (Exchange Calculation), MD run file preparation from ```gmx grompp``` (Simulation Preparation), the actual MD simulation (MD Simulation), and the sum of the previous three (Timestep Total). Another larger plot is created at right showing the cumulative time of the simulation in hours. In general simESI-fmm is fairly well optimized. While performance is dependent on both CPU/GPU strength, in general, the calculations done by simESI-fmm take <5% of total compute time. The majority of computational expense is from the actual MD simulation or ```gmx``` calls including ```grompp``` or ```pdb2gmx```.
 
<p align="center">Example output from runtime.py using ubiquitin as an example.</p>

![runtime.py Example](/assets/runtime_ex_ubq.png) 
