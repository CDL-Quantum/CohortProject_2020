# Ising data set
The provided data set is a collection of samples from an Ising Hamiltonian with pairwise interactions.

# Files

The following two files are provided:

* Training set: `train_{hashcode}.csv`
* Validation set: `val_{hashcose}.csv`
* List of pair of correlated variables: `correlated_features.csv`

You may use the first file to train a Boltzmann machine and use the second data set to validate the generative model.

# Tasks

1. Generate 10000 samples
2. Ensure the samples are generated from an equilibrated model. What metric do you use to prove equilibrium?
3. Report the mean of each Ising variable
4. Report the mean correlation of pairs of variables listed in `correlated_features.csv`