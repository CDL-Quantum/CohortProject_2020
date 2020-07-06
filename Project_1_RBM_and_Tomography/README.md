![CDL 2020 Cohort Project](../figures/CDL_logo.jpg)
## Project 1: Machine Learning

This project will guide you through using a machine learning algorithm -- the Restricted Boltzmann machine (RBM) -- on molecular data. If you want some more machine learning exposure before you begin your tasks, supplemental reading on RBMs can be found [here](\href{https://qucumber.readthedocs.io/en/stable/_static/RBM_tutorial.pdf}), and a guide on training RBMs on a dummy dataset has been provided in a [Jupyter notebook](https://github.com/CDL-Quantum/CohortProject_2020/blob/master/Project_1_RBM_and_Tomography/RBM_train_dummy_dataset.ipynb). The base RBM code is located in [RBM_helper.py](https://github.com/CDL-Quantum/CohortProject_2020/blob/master/Project_1_RBM_and_Tomography/RBM_helper.py). We've intentionally left it to be bare-bones and easy to use.

Open up [Project1_LandingPage.pdf](https://github.com/CDL-Quantum/CohortProject_2020/blob/master/Project_1_RBM_and_Tomography/Project1_LandingPage.pdf) to begin learning about your tasks for this week!

## Tasks include:
* Train an RBM to reconstruct the potential energy of moleculary Hydrogen as the distance between the H atoms changes.
* Train an RBM to reconstruct the groundstate of a Rydberg atom quantum computer.

## Further Challenges:
* Additional data has been provided for 
[LiH and BeH2.(https://github.com/CDL-Quantum/CohortProject_2020/tree/master/datasets/qubit_molecules).
Consider an RBM reconstruction for these more complicated molecules.  What is the fundamental difference between these and the systems studied in the above Tasks?
* Explore other unsupervised machine learning techniques besides RBMs, such as recurrent neural networks (RNNs). Which can be used to represent molecular wavefunctions?
* What other quantum wavefunction datasets can you obtain? Is a “standard dataset” for machine learning quantum chemistry possible?
* D-Wave has included some data on a 
[hard optimization problem](https://github.com/CDL-Quantum/CohortProject_2020/tree/master/datasets/IsingSamplesDW). 
Can you use your machine learning skills to uncover the correlations between the Ising variables?


## Business Application
For each week, your team is asked to complete a Business Application. Questions you will be asked are:

* Explain to a layperson the technical problem you solved in this exercise.
* Explain or provide examples of the types of real-world problems this solution can solve.
* Identify at least one potential customer for this solution - ie: a business who has this problem and would consider paying to have this problem solved.
* Prepare a 90 second video explaining the value proposition of your innovation to this potential customer in non-technical language.

For more details refer to the [Business Application found here](./Business_Application.md)
