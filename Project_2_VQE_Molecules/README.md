![CDL 2020 Cohort Project](../figures/CDL_logo.jpg)
## Project 2: VQE
Within this module, you will explore the Variational Quantum Eigensolver (VQE) for
constructing potential energy surfaces for small molecules.
A very brief introduction to the main ideas behind the VQE are 
[here.](https://github.com/CDL-Quantum/CohortProject_2020/blob/master/CDL_2020_docs.pdf)
Open up [Project2_LandingPage.pdf](https://github.com/CDL-Quantum/CohortProject_2020/blob/master/Project_2_VQE_Molecules/Project_2_LandingPage.pdf)
to begin learning about your tasks for this week!

## Tasks and Challenges
Using H2 and H2O as examples, you will implement the main steps of the VQE process from the setup of the Hamiltonian to obtaining the circuit for the IBM quantum computer.  A number of additional challenges and 
business-related questions are suggested in the above.

We have implemented the main steps of the VQE process in the following notebooks :

[Classical Method Solution](https://github.com/tina-seb/CohortProject_2020/blob/master/Project_2_VQE_Molecules/S1_Classical_Methods_Demo.ipynb)

[Hamiltonian Generation Summary](https://github.com/tina-seb/CohortProject_2020/blob/master/Project_2_VQE_Molecules/S2_Hamiltonian_gen_Summary.ipynb)

After running the original h2 example we go over the molecules='h2', 'h4', 'lih', 'h2o', 'n2', 'nh3', basies='sto-3g','6-31g', and qubit_transfms='jw','bk' to assess the scope of the final measurements  by printing out information on the QWC fragments and do Random test of the Unitary ... we suppress printing information for more than 9 mutually commuting fragments.

[Unitary Ansatz](https://github.com/tina-seb/CohortProject_2020/blob/master/Project_2_VQE_Molecules/S3_Unitary_Ansatz_H2.ipynb)

[Measurement Summary](https://github.com/tina-seb/CohortProject_2020/blob/master/Project_2_VQE_Molecules/S4_Measurement_Summary.ipynb)

[Circuits](https://github.com/tina-seb/CohortProject_2020/blob/master/Project_2_VQE_Molecules/S5_Circuits-H2_on_IBMq-sussex.ipynb)

## Business Application
Once again, your team is asked to complete a Business Application. Questions you will be asked are:

* Explain to a layperson the technical problem you solved in this exercise.
* Explain or provide examples of the types of real-world problems this solution can solve.
* Identify at least one potential customer for this solution - ie: a business who has this problem and would consider paying to have this problem solved.
* Prepare a 90 second video explaining the value proposition of your innovation to this potential customer in non-technical language.

For more details refer to the [Business Application found here](./Business_Application.md)

## Presenting your results in your pull request
For your pull request, consider the following for the presentation of your final results:
- Work entirely in the directory for Project 2.
- Edit this README.md file with a highlight of your main technical results.  Provide links to any other files with your detailed results, e.g. Jupyter notebooks.
- For your Business Application, feel free to provide your answers directly in the 
[Business_Application.md](./Business_Application.md) file.
- Do not directly upload your video file (or any other large files) to the repository.  Instead, provide a link e.g. to a YouTube video, or a Google Drive file.
- Include a file contributions.md that lists the contributions of each group member.
