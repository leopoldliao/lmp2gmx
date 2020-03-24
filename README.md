#  LAMMPS data file to GROMACS `.gro` & `.itp` forcefield translator

User must specify input and output file names in `in.lmp2gmx` file. 

## Authors 
* **[Chris Lorenz](chris.lorenz@kcl.ac.uk)**, 
* **[Mohamed Ali al-Badri](mohamed.al-badri@kcl.ac.uk)**

This script assumes the following: 

* A very particular format of lammps data file -- use `data_to_data.py` in Robbie's GOAzide repo.
* A list of atom names in the `Masses` section of the LAMMPS data file -- these are output to go in the [ atomtypes ] section of the gromacs topology.
* The [ pairs ] section of the gromacs topology are output separately -- insert this in the main `.itp` file above the [ dihedrals ] section.

