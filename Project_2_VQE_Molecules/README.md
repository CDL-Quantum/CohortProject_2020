![CDL 2020 Cohort Project](../figures/CDL_logo.jpg)
## Project 2: VQE
Within this module, you will explore the Variational Quantum Eigensolver (VQE) for
constructing potential energy surfaces for small molecules.
A very brief introduction to the main ideas behind the VQE are 
[here.](https://github.com/CDL-Quantum/CohortProject_2020/blob/master/CDL_2020_docs.pdf)
Open up [Project2_LandingPage.pdf](https://github.com/CDL-Quantum/CohortProject_2020/blob/master/Project_2_VQE_Molecules/Project_2_LandingPage.pdf)
to learn in detail about the project.

## Tasks and Challenges
We study classical as well as quantum methods to construct Potential Energy Surface (PES) for several small molecules like H2, H20, LiH, N2. For classical simulation we use standard methods like CCSD, HF, FCI and for quantum simulation we use VQE. The work is presented in following notebooks as follows:

[Classical Method Solution](https://github.com/tina-seb/CohortProject_2020/blob/master/Project_2_VQE_Molecules/S1_Classical_Methods_Demo.ipynb)

In this notebook we obtain the electronic PES in minimal basis (STO-3G) for different molecules using classical methods like FCI, HF, CISD, CCSD available in the tequila package. We show PES for several molecules -H2, H20, LiH, N2 obtained by fore-mentioned methods.

[Hamiltonian Generation Summary](https://github.com/tina-seb/CohortProject_2020/blob/master/Project_2_VQE_Molecules/S2_Hamiltonian_gen_Summary.ipynb)

After running the original H2 example we try Hamiltonian generation for other molecules 'H4', 'LiH', 'H2O', 'N2', 'NHh3'in both 'sto-3g','6-31g' basis under Bravyi-Kitaev and Jordan Wigner Transformations to assess the scope of the final measurements  by printing out information on the Qubit Wise Commuting fragments and do Random test of the Unitary.

[Unitary Ansatz](https://github.com/tina-seb/CohortProject_2020/blob/master/Project_2_VQE_Molecules/S3_Unitary_Ansatz_H2.ipynb)

[Measurement Summary](https://github.com/tina-seb/CohortProject_2020/blob/master/Project_2_VQE_Molecules/S4_Measurement_Summary.ipynb)

[Circuits](https://github.com/tina-seb/CohortProject_2020/blob/master/Project_2_VQE_Molecules/S5_Circuits-H2_on_IBMq-sussex.ipynb)

In this notebook we run, using the data obtained from previous notebooks, the experiment on a 5 qubit quantum computer (Sussex) from IBMq. The results show that noise affect considerably. However, it manages to map a considerably good curve and as soon as quantum hardware gets better it is expected to get better results.

## Business Application
For the Business Application we consider several plausible applications. For instance:
- Water Treatment 

- Fertlizer Manufacturing 

- Efficient Battery Production 

- Small Rocket Propellants 

For more details refer to the [Business Application found here](./Business_Application.md)

We also made a [video]() in which we explain what can be expected of all these applications nowadays and in the near future.

## Future work and Possibilities

We intend to explore and implement VQE for obtaining excited states and use VQE subroutines to obtain additional properties like charge density, dipole moment etc. on various hardwares and platforms available today. Additionally, we also consider possibility of implementing and benchmarking ADAPT-VQE[10] for exact simulations on current hardware. 

**//Delete this subheading below and info under it once done editing readme//**
## Presenting your results in your pull request
For your pull request, consider the following for the presentation of your final results:
- Work entirely in the directory for Project 2.
- Edit this README.md file with a highlight of your main technical results.  Provide links to any other files with your detailed results, e.g. Jupyter notebooks.
- For your Business Application, feel free to provide your answers directly in the 
[Business_Application.md](./Business_Application.md) file.
- Do not directly upload your video file (or any other large files) to the repository.  Instead, provide a link e.g. to a YouTube video, or a Google Drive file.
- Include a file contributions.md that lists the contributions of each group member.
