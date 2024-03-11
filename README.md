# Pathogen Sensing Project

This is the companion repository for the paper "Shedding light on bacterial fitness in a tug-of-war with liquid crystal emulsions" containing simulation, analysis, and visualization code.

<img src="https://github.com/bzjan/pathogen_sensing/assets/5402654/6097459d-43b7-4d79-83db-cbbb201d32f7" width="400">


### 1. System Requirements
#### Hardware
* OS: Windows 11 10.0.22631
* Processor: Intel i7-13700F
* RAM: 64 GB
* Storage: ~ 1 GB

#### Software
* [Mathematica](https://www.wolfram.com/mathematica/) (or free [Wolfram Engine](https://www.wolfram.com/engine/)) 13.3
* [Matlab](https://www.mathworks.com/products/new_products/latest_features.html) 2023b
* gcc 13.2.0 via [msys2](https://www.msys2.org/)
* [Visual Studio Code](https://code.visualstudio.com/download) 1.87
* [CMake Tools extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode.cmake-tools) 1.17.17

### 2. Installation guide
#### Get this repository with e.g.
```
git clone https://github.com/bzjan/pathogen_sensing.git
```

#### Install the C++ simulator:
The following instructions will take about 5 minutes:
* For a free C++ compiler, download and install msys2 from: https://www.msys2.org/
* Open msys terminal via: `windows key > msys` and enter:
```
# update packages
pacman -Syu
```
After installation agree to close the msys terminal and start a new one via: `windows key > msys` and enter:
```
# update packages
pacman -Su

# gcc, make, cmake, boost, mpi, gdb, eigen
pacman -S mingw-w64-x86_64-gcc mingw-w64-x86_64-make mingw-w64-x86_64-cmake mingw-w64-x86_64-boost mingw-w64-x86_64-msmpi mingw-w64-x86_64-gdb  mingw-w64-x86_64-python-pygments mingw-w64-x86_64-eigen3
```
* In Visual Studio Code install and use the CMake Tools extension to automatically detect the installation files, build the simulator executable (takes about 5s) and run unit tests.
* The compiled executable must be placed in a subfolder `build` in the folder where the Mathematica files reside for Mathematica to be able to call the executable.

### 3. Demo
* To test the C++ simulation code, we included unit tests. They can be executed with the CMake Tools extension in Visual Studio Code. When run, all unit tests pass successfully.
* To test the simulation, you need to have built the simulator executable. Afterwards open `Mathematica/run_simulations.wl` in Mathematica and go to subsection `3d position with 3d orientation vector/scan: single simulation`. To initialize the required functions use the `Run package` button in the top left corner. Execute the code block with `Shift+Enter` to start the simulation. It should take about 10s.
The output will be available in Mathematica for visualization and analysis, which can be done in the following block (10s). The resulting plot will be located in output/test
* Expected run times for MATLAB code:
   * dropletLightProp.m: Approximately 1-2s (with program options set to false)
   * trackDroplets.m: Approximately 10 minutes to process 9000 frames containing about 30 droplets,
   * frictionCalibrationForTrackedDroplets.m: Approximately 50s to process 22 droplets
* Datasets to demo MATLAB code:
   * dropletLightProp.m: test_dropletTilt_thetaPhi_analyticDirector_wvl_550_DnLC_0.18, located within Data folder
      * Modify variables thetavals and phivals (located on lines 16 and 17 respectively) to both be 0
      * Set showImagesLiveQ (line 33) to true, and ensure that other program options (lines 32, 34-38) are set to false
      * This dataset contains a single electromagnetic field generated from Jones calculus for an upright droplet (theta = 0, phi = 0)
   * trackDroplets.m: Load Supplementary Video 4 to extract droplet tilts
   * frictionCalibrationForTrackedDroplets.m: Use dataset extracted by trackDroplets.m to obtain friction coefficient

### 4. Instructions for use
The repository contains files for 
* analyzing experimental data
  * Mathematica/droplet_tracker.wl (extract droplet motion and light intensities over time)
  * Mathematica/droplet_optics.wl (Computer liquid crystal arrangement and simulate light transmission through droplets using Jones calculus)
  * Matlab/frictionCalibrationForTrackedDroplets.m (measure effective droplet friction coefficient)
  * Matlab/dropletLightProp.m (Propagation of electromagnetic waves beyond the droplet)
  * Matlab/trackDroplets.m (Measuring droplet tilt over time)
* orchestrating and analyzing simulations
  * Mathematica/run_simulations.wl (define parameters for simulation scans, analyze and visualize results)
* creating the manuscript figures
  * Mathematica/figures.wl

To recreate the manuscript figures you need to open figures.wl in Mathematica. Adapt the paths in subsection `Constants/parametersPaths: paths` and execute the notebook with the button `Run package`. Creating the figures will take about 10 minutes.
