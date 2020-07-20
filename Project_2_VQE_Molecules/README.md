![CDL 2020 Cohort Project](../figures/CDL_logo.jpg)
## Project 2: VQE

## Tasks and Challenges
We study classical as well as quantum methods to construct Potential Energy Surface (PES) for several small molecules like H2, H20, LiH, N2. For classical simulation we use standard methods like CCSD, HF, FCI and for quantum simulation we use VQE. The work is presented in following notebooks as follows:

[Classical Method Solution](https://github.com/tina-seb/CohortProject_2020/blob/master/Project_2_VQE_Molecules/S1_Classical_Methods_Demo.ipynb)

In this notebook we obtain the electronic PES in minimal basis (STO-3G) for different molecules using classical methods like FCI, HF, CISD, CCSD available in the tequila package. We show PES for several molecules -H2, H20, LiH, N2 obtained by fore-mentioned methods.

[Hamiltonian Generation Summary](https://github.com/tina-seb/CohortProject_2020/blob/master/Project_2_VQE_Molecules/S2_Hamiltonian_gen_Summary.ipynb)

After running the original H2 example we try Hamiltonian generation for other molecules 'H4', 'LiH', 'H2O', 'N2', 'NHh3'in both 'sto-3g','6-31g' basis under Bravyi-Kitaev and Jordan Wigner Transformations to assess the scope of the final measurements  by printing out information on the Qubit Wise Commuting fragments and do Random test of the Unitary.

[Unitary Ansatz](https://github.com/tina-seb/CohortProject_2020/blob/master/Project_2_VQE_Molecules/S3_Unitary_Ansatz_H2.ipynb)

Running the and checking the Unitary ansatz for H2 in two possible bases.
We identified problem with the energy at R=1.5A the two bases give different energies!

[Measurement Summary](https://github.com/tina-seb/CohortProject_2020/blob/master/Project_2_VQE_Molecules/S4_Measurement_Summary.ipynb)

After running the original H2 example we go over the molecules: 'H2', 'H4', 'LiH', 'H2O', 'N2', 'NH3', the two bases 'sto-3g' and '6-31g', and the qubit_transfms 'jw' and 'bk' to assess the scope of the final measurements by printing out information on the QWC fragments and doing random tests of the Unitary. 
We suppress the printing of information for more than 9 mutually commuting fragments.

[Circuits](https://github.com/tina-seb/CohortProject_2020/blob/master/Project_2_VQE_Molecules/S5_Circuits-H2_on_IBMq-sussex.ipynb)

In this notebook we run, using the data obtained from previous notebooks, the experiment on a 5 qubit quantum computer (Sussex) from IBMq. The results show that noise affect considerably. However, it manages to map a considerably good curve and as soon as quantum hardware gets better it is expected to get better results.

[Research Files Readme](https://github.com/tina-seb/CohortProject_2020/blob/master/Project_2_VQE_Molecules/Research_Files/README.md)

We elaborate in more detail the results of our calculations and the problems we encountered in the process.

## Business Application
For the Business Application we consider several plausible applications. For instance:
- Water Treatment 

- Fertlizer Manufacturing 

- Efficient Battery Production 

- Small Rocket Propellants 

For more details refer to the [Business Application found here](./Business_Application.md)

We also made a [video](https://www.youtube.com/watch?v=__A3Da354DE&feature=youtu.be) in which we explain what can be expected of all these applications nowadays and in the near future.

## Future work and Possibilities
The first task is to implement noise mitigation of the q-devices, using their calibration, and to test the range variation of the results.
[This is critical for VQE use.](https://github.com/tina-seb/CohortProject_2020/blob/master/Project_2_VQE_Molecules/Research_Files/README.md)
We also intend to explore and implement VQE for obtaining excited states and use VQE subroutines to obtain additional properties like charge density, dipole moment etc. on various hardware and platforms available today. Additionally, we also consider the possibility of implementing and benchmarking ADAPT-VQE[10] for exact simulations on current hardware. 
