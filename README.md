# vsc-heom

A HEOM code for computing the rate constants of a model molecular system in contact with an optical cavity as described in [Quantum Dynamics of Vibrational Polariton Chemistry]{https://arxiv.org/abs/2210.05550}.  Here a number of different programs that treat different configurations of the system and using different initial conditions for evaluating the rate constant.  Currently the CMake build script is setup to compile the five source files (.cpp into five separate executables) for evaluating using the symmetrized side operator for evaluating the side-side correlation function:
***
* rate_theory/main.cpp -> heom_mol.x: A reaction coordinate coupled to an unstructured Debye bath. 
* rate_theory/main_cavity.cpp heom ->_cavity.x: A reaction coordinate coupled to an unstructured Debye bath and to a perfect (no loss) cavity
* rate_theory/main_cavity_loss.cpp -> heom_cavity_loss.x: A reaction coordinate coupled to an unstructured Debye bath and to a lossy cavity (modelled by attaching a Debye bath to the cavity degrees of freedom)
* rate_theory/main_mol_2d.cpp -> heom_mol_2d.x: A reaction coordinate coupled to a structured bath (Debye + underdamped mode)
* rate_theory/main_cavity_loss_2d.cpp -> heom_cavity_loss_2d.x: A reaction coordinate coupled to a strcutured bath (Debye + underdamped mode) and to a lossy cavity

Additional source files are present that compute an alternative correlation function obtained by first performing a thermalisation run of the HEOM these are present in:
***
* rate_theory/main_thermalise.cpp -> heom_mol_thermalise.x: A reaction coordinate coupled to an unstructured Debye bath. 
* rate_theory/main_cavity_loss_thermalise.cpp -> heom_mol_thermalise.x: A reaction coordinate coupled to an unstructured Debye bath. 

Finally source files are included for evaluating the IR spectrum in a number of different cases:
***
* ir_spectrum/main_ir.cpp -> heom_mol_ir.x: A reaction coordinate coupled to an unstructured Debye bath. 
* ir_spectrum/main_cavity_loss_ir.cpp -> heom_cavity_loss_ir.x : A reaction coordinate coupled to an unstructured Debye bath and to a lossy cavity (modelled by attaching a Debye bath to the cavity degrees of freedom)
* ir_spectrum/main_mol_ir_2d.cpp -> heom_mol_2d_ir.x: A reaction coordinate coupled to a structured bath (Debye + underdamped mode)
* ir_spectrum/main_cavity_loss_ir_2d.cpp -> heom_cavity_loss_2d_ir.x: A reaction coordinate coupled to a strcutured bath (Debye + underdamped mode) and to a lossy cavity

Currently the thermalisation code and IR spectrum code will not be compiled, however, by setting the cmake variables BUILD_THERMALISE and BUILD_IR to On these files will be added to the build process.

## Disclaimer
This software is provided to support the findings presented in [Quantum Dynamics of Vibrational Polariton Chemistry]{https://arxiv.org/abs/2210.05550} and is very much a research code.  As such, it is currently very much specialised to perform the calculations presented in this work.  While the underlying HEOM method is considerably more general than the application considered here, and even for this application other observables could readily be evaluated, no effort has been made to expose this functionality through the input files.  In principle, this code can be modified to treat arbitrary system bath models coupled to Debye baths, but at the current moment this is not something that will be actively pursued with this code base. 

## Dependencies
External:
- Required: [RapidJSON](https://rapidjson.org/) input file parsing


Submodules (automatically downloaded if the respository is cloned using git clone --recurse-submodules https://github.com/llindoy/vsc-heom.git):
- [linalg](https://github.com/llindoy/linalg) wrapper for the BLAS and LAPACK libraries

# Compile Instructions
This code requires cmake version 3.11 in order to compile. From hyperfine moment base directory (${vsc-heom_base}) run:
```console
mkdir build
cd build
cmake ../
make
make install
```

This builds the executable as ${vsc-heom_base}/bin/${executable_name}, where the executable names are outlined above.

To build the executables that use a thermalised initial condition modify the cmake command to:
```console
cmake -DBUILD_THERMALISE=On ../
```

To build the executables for evaluating the IR spectrum modify the cmake command to:
```console
cmake -DBUILD_IR=On ../
```

Both options can be specified simultaneously.

The cmake build command will scan the directory (${vsc-heom_base}/external) for any of the required external libraries, and use them if present.  If they aren't it will attempt to download the github repositories and make them available using the cmake FetchContent package.

## Running the Software
Each of the generated executables expects a string specifying the path to the input file and will write the output to stdout and some additional information to stderr.  A series of example input files and the resultant output files are provided in the examples directory, further information on the parameters in these input files are provided in the input_file_structu   re.in file in each of the example subdirectories.  Further scripts for evaluating the rate constants and IR spectra from the output correlation functions are provided as extract_rates.py and plot_spectrum.py. 

The input files for the "thermalise" variants of each approach only differ from the standard input files in that it is necessary to include an additional argument "nequil" = ${integer number of equilibration steps}.  If this argument is included in the input files used for the standard code it is ignored.
