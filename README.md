# vsc-heom

A HEOM code for computing the rate constants of a model molecular system in contact with an optical.  Here a number of different programs that treat different configurations of the system and using different initial conditions for evaluating the rate constant.  Currently the CMake build script is setup to compile the five source files (.cpp into five separate executables) for evaluating using the symmetrized side operator for evaluating the side-side correlation function:
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

Currently the thermalisation code and IR spectrum code will not be compiled, however, by uncommenting the relevant lines in the src/CMakeLists.txt file it is possible to add these files to the build process. 


## Dependencies
External:
- Required: [RapidJSON](https://rapidjson.org/) input file parsing


Submodules (automatically downloaded if the respository is cloned using git clone --recurse-submodules https://github.com/llindoy/vsc-heom.git):
- [linalg](https://github.com/llindoy/linalg) wrapper for the BLAS and LAPACK libraries

# Compile Instructions
This code requires cmake version 3.11 in order to compile. From hyperfine moment base directory (${hyperfine_base}) run:
```console
mkdir build
cd build
cmake ../
make
make install
```

This builds the executable as ${vsc-heom_base}/bin/${executable_name}

where the executable names are outlined above.

The cmake build command will scan the directory (${vsc-heom_base}/external) for any of the required external libraries, and use them if present.  If they aren't it will attempt to download the github repositories and make them available using the cmake FetchContent package.

## Running the Software
Each of the generated executables expects a string specifying the path to the input file and will write the output to stdout and some additional information to stderr.  A series of example input files and the resultant output files are provided in the examples directory, further information on the parameters in these input files are provided in the input_file_structure.in file in each of the example subdirectories.  Further scripts for evaluating the rate constants and IR spectra from the output correlation functions are provided as extract_rates.py and compute_spectrum.py.
