# Method-development2(Kirk-Wood buff theory)

This repository consist Python codes written for calculating preferentail intercation coefficients using Kirk-wood-Burf integral apparoch. This methods are useful to calculate preference of solvent towrads the solute. example of urea prefernce towards tryptophan moleucles is shown.


## Getting Started

These instructions will get you a copy of the project up and running on your local machine for testing purposes. 

### Prerequisites

Python libraries required to use the code

```
sys
math
matplotlib
scipy
numpy 
savitzy_golay 
copy

```

### Installing
Python 3 or above version is required to setup the dependencies. This code uses Savitzky_golay filter

```
use commands below to activate loacal environment anaconda
conda create -n method python=3.6
source activate method

e.g. to install python libraries listed above

conda install -c omnia matplotlib  (similarly install other libraries using conda or pip)



## Running the codes

KB_integral_processing.py takes 3 standard inputs from user i.e. radial distribution function files for solvents in your solutaion (obtained from the trajectory file processing in CHRAMM simulation program), pdb file, model residue name 

To run the automated tests for this system
kb.sh can be used as:
bash kb.sh



