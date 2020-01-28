# DiskDyn

An Nvidia/CUDA C based simulator for the dynamical evolution of debris disks.

DiskDyn can models random swarms of gravitating bodies in a disk or cloud (N-body code like). It also produces realistic images of the systems and calculates their full SEDs, including scattered and thermal light, taking into account the scattering angles. 

# Installing

The github package contains the source code, test files, tools to calculate optical tables, and a cfitsio directory. The tools to calculate optical tables are required to pre-calculate optical constants for DiskDyn. Pre-calculated tables for astronomical silicates, iron, olivine, organics, orthopyr, troilite, and water ice can be downloaded from the page website at http://merope.as.arizona.edu/~agaspar/d3d.html

Once downloaded, cfitsio needs to be compiled first, unless the system already has cfitsio installed.

$ cd cfitsio
$ ./configure
$ make
$ cd ../

If a system-wide cfitsio is used, the Makefile in src needs to be edited to point to proper libraries and include locations. Otherwise, a simple make command will compile all the necessary components.

$ make

# Running DiskDyn

DiskDyn takes an input parameter file that contains all the information it needs to run or a previous data output file. Example parameter files can be found in the tests directory.

$ ./bin/ddyn -i tests/SolarSystem.param

The variables are explained in the parameter file.

The optical constants need to be placed somewhere and access path to them given in the parameter file. The example parameter files have them living in a directory called "optical".

All file locations are given relative to the DiskDyn directory. This may present issues for servers where file paths are given relative to the home directory. Contact me if that's an issue for you. A correction for this will be given soon.

# Calculating new optical constants for DiskDyn

If new optical constants are required, they need to be calculated using the runmie program located in the tools directory. The tool requires l,n,k values. The number of dust sizes, scattering angles, and the minimum and maximum dust sizes can be given to the tool, if the default values are not adequate.

# Example simulations using DiskDyn

DiskDyn calculates the dynamical evolution of debris disks and dust clouds, dynamically perturbed by planets and/or gravitating swarm of bodies. The code also calculates images in both scattered light and thermal emission at requested wavelengths (assuming mie scattering) accounting for scattering angles. Additionally, complete SEDs are also calculated at requested time intervals.
