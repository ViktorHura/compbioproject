
# Genetic Algorithm approach for finding hidden motifs in DNA sequences

This is the source code for our project for `1000WETCOB Inleiding tot computationele biologie` 


## Installation

Install this project by installing the packages in the requirements file to your python environment.

```bash
  pip install -r requirements.txt
```

Note that we used the `colorama` package for visualisations in the console, which we have only tested on Windows.
    
## Structure

`data/` contains the DNA sequence files that we used for testing the algorithms, 
along with a simple script to extract the UCSC-Cat sequences from the provided fasta file.

`output/` is where our scripts will write output csv files.

`src/` contains the source code for our genetic algorithm, we wrote it in a modular way for flexibility.

`util.py` contains functions for generating PSSM and calculating PSSM score.

`gibbs.py` is the script for the Gibbs motif sampling algorithm.

`genetic.py` is the script for our genetic algorithm approach.


## Running Tests

To run tests, see the first lines of the `gibbs.py` and `genetic.py` scripts.
They contain configuration variables which you can change to your needs.

