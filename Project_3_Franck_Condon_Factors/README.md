![CDL 2020 Cohort Project](../figures/CDL_logo.jpg)
## Project 3: Franck-Condon Factors

In this project we calculate Franck-Condon Factors, which is related to the intensities in spectra of vibrational transitions across electronic surfaces. Spectra can be measured through experiments, but having accurate theoretical calculations of these Franck-Condon Factors allow scientists to instead predict these experimental results. This is especially useful if the chemical species is expensive, difficult to acquire or difficult to study or if the experiment itself is expensive or difficult.

A very brief introduction to the main ideas behind the project are
[here.](https://github.com/CDL-Quantum/CohortProject_2020/blob/master/CDL_2020_docs.pdf)

In the [Project3_LandingPage.pdf](https://github.com/CDL-Quantum/CohortProject_2020/blob/master/Project_3_Franck_Condon_Factors/Project3_LandingPage.pdf), a more technical information about the chemistry and physics behind these Franck-Condon Factors is given. There, you can also find the descriptions of whose solutions are provided below.

# Tasks
This whole project is about theoretical quantum chemistry and spectroscopy in particular.  
On the technical side, we are studying different theoretical methods to investigate and predict chemical properties of molecules. In particular we are focusing on the so called the vibronic transitions and Franck-Condon factors of molecules.
**Task 1:**   
In this task we use a very simple theoretical model to study the properties of the simplest molecule in the world, **molecular Hydrogen**. The purpose of this task is twofold - a) this is educational, and any newcommer to the field can quickly get familiar with the technical basics; b) In this task we show, that the basic theoretical model (harmonic oscillator approximation) used throughout the project is sound. To this end, we do the theoretical calculations for hydrogen, then compare the results with actual experimental data, and verify that the theory gives correct predictions.  
The calculations for first task for the Hydrogen molecule can be found in [Task1 jupyter notebook](https://github.com/hay-k/CohortProject_2020_w3g7/blob/master/Project_3_Franck_Condon_Factors/Task1.ipynb). Particularly, we calculate transitions from n=0 state of H<sub>2</sub> to 10 other vibronic levels for H<sub>2</sub><sup>+</sup>. The figures below show the simulation results from our code implemented in `FCF_helper.py` file which does all the calculations, and the snapshot below is from the paper by Berkowitz and Spohr, _Journal of Electron Spectroscopy and Related Phenomena_, **2**(2):143–152 (1973). The points or vertical bars in our plots correspond to the peak values in the original figure. We can see the great resemblance of our result with the original paper.                                                                              
![](Vibronic_spectrum.png "our catterplot") ![](Vibronic_sticks.png "our plot: sticks")  
<img src="Berkowitz.png" width="850"/> 

**Task 2:**  
In this task we introduce a more sophisticated technical tool (FC.cxx), which can calculate the same chemical properties of molecules as in Task 1, but it can do it for a wide range of molecules, not only hydrogen. This tool is created by a huge research effort in a university, and has been extensively tested already and verified with experimental results. We do not need to test the tool. We just use this tool to investigate the properties of the **Vanadium3 (V<sub>3</sub>) molecule**. For the visualizations of the results output from FC.cxx tool see our [Task2 jupyter notebook](https://github.com/hay-k/CohortProject_2020_w3g7/blob/master/Project_3_Franck_Condon_Factors/Task2.ipynb).  

**Task 3:**  
In this task we introduce yet another sophisticated technical tool Gaussia Boson Sampler (GBS) - a photonic special-purpose sampling device by Xanadu programmed via their Strawberryfields software ramework. GBS which tackles the problem from a completely different viewpoint than the tool in the previous task. We investigate the properties of the V3 molecule with this tool and verify that the results are in line with the results produced by the tool in the previous task. One important advantage of this new tool compared to the previous one is that it runs significantly faster for large molecules, and enables theoretical investigation of a few molecules, which are very time consuming and difficult with the previous tool. The investigations are visualized in our [Task3 jupyter notebook](https://github.com/hay-k/CohortProject_2020_w3g7/blob/master/Project_3_Franck_Condon_Factors/Task3.ipynb).   

**Challenge 1:** 
In this challenge we further investigate the tool introduced in Task3 to reveal more cons and pros about it. This time calculations are done using the loop hafnian approach. The loop hafnian approach uses GBS which would allow the Franck-Condon factors to be calculated using a quantum circuit. We use the result of Task 3 to provide data to a skeleton code provided that uses loop hafnians to calculate the Franck-Condon Factors.

**Challenge 2:** 
In this challenge we are given a freedom to pick a molecule of our liking, and investigate its properties. As discussed during the call, this is a good point for business proposal entrance, so we need to choose a nice molecule for investigation, which can have some interesting business applications

**Challenge 3:** 
We analyze and report all the advantages and disadvantages of all tools used above.

**Challenge 4:** 
We investigate the advantages and disadvantages of codes licensed for the public domain and those that are licensed for private use

# Business Application
One again, your team is asked to complete a Business Application. Questions you will be asked are:

* Explain to a layperson the technical problem you solved in this exercise.
* Explain or provide examples of the types of real-world problems this solution can solve.
* Identify at least one potential customer for this solution - ie: a business who has this problem and would consider paying to have this problem solved.
* Prepare a 90 second video explaining the value proposition of your innovation to this potential customer in non-technical language.

For more details refer to the [Business Application found here](./Business_Application.md)

## Presenting your results in your pull request
For your pull request, consider the following for the presentation of your final results:
- Work entirely in the directory for Project 3.
- Edit this README.md file with a highlight of your main technical results.  Provide links to any other files with your detailed results, e.g. Jupyter notebooks.
- For your Business Application, feel free to provide your answers directly in the 
[Business_Application.md](./Business_Application.md) file.
- Do not directly upload your video file (or any other large files) to the repository.  Instead, provide a link e.g. to a YouTube video, or a Google Drive file.
- Include a file contributions.md that lists the contributions of each group member.
