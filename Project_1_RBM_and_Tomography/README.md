# CDL-2019-1st-Project

## Project Restricted Boltzmann Machines (RBM) as a state ansatz (Patrick Huembeli)

### Additional reading to the Project

- There are many introdcutions into RBMs [this](https://arxiv.org/abs/1806.07066) is a quite extensive one.
- For more explanations about the derivation of the update rule of the RBM go [here](https://qucumber.readthedocs.io/en/stable/_static/RBM_tutorial.pdf)
- For more info about [Hopfield networks](https://page.mi.fu-berlin.de/rojas/neural/chapter/K13.pdf), the XOR problem (page 347), the Hebbian learning and why the data becomes minimum energy (page 354) have a look [here](https://page.mi.fu-berlin.de/rojas/neural/chapter/K13.pdf).
- To understand the connection between quantum states and RBMs have a look at [this](https://arxiv.org/abs/1703.05334) and [this](https://arxiv.org/abs/1606.02318) paper.

## Get started

You can get started the following way:
- The 'RBM' folder contains a fully working restricted Boltzmann machine and a dummy dataset that only contains strings of [1,0,1,0,...] and [0,1,0,1,...] which can be used as a training set for first tests. If you run `RBM_train_dummy_dataset.ipynb` it should train on this simple dummy set.
- Make the file `RBM_train_dummy_dataset.ipynb` work on your computer and try to understand the code.
- The whole RBM is defined in `RBM_helper.py` you will have to understand this file.
- Test if the training worked well by sampling from the trained RBM. If you sample after training on this dummy data you should only get the strings [1,0,1,...] and [0,1,0,...] because this is the only data the RBM has seen. If you get other strings, your RBM is not trained perfectly. This will happen at some point and this is totally normal. Try to improve the training and play with this file.

**First assignment:**
- Understand the code in `RBM_helper` as good as possible and write a pseudocode of it.
- Examples for pseudo codes can be found e.g. [here](https://en.wikibooks.org/wiki/LaTeX/Algorithms) or [here](https://tex.stackexchange.com/questions/163768/write-pseudo-code-in-latex)
- You can do it by hand as a sketch.

**Second assignment a):**
- Train the RBM on the Bars and Stripes OR MNIST data.
- To start I recommend Bars and Stripes, because it needs much less computational resources. In the file `Bars_and_Stripes.py` you can see how the data set is generated and saved to a `.npy` file.
- The bars and stripes "images" are only 4x4. Later you should scale them to bigger images.
- We have to convert the 4x4 pixel images into a 16 dimensional vector and put them into a RBM with the same input dimension. Find out how to do this with for example 'numpy'.
- Train the RBM the same way you did it on the dummy data set from before just with a different input dimension.
- Sample from it and see what you obtain. Document this in the final report. Describe what happens. Are the bars and stripes you sample good or are they blurry or just random noise? What happens if you increase the training epochs or the steps of the Gibbs sampling? Be aware that you have to transform the output of the RBM back to a 4x4 image. The same way you transformed the 4x4 image before to a 16 dimensional vector.
- Explore what happens if you change the number of hidden units. In the code so far we set the number of visible and hidden units equal. Can we still learn anything if we take e.g. 4 times less hidden units? How do the images look like? Explore and discuss! How does this scale with the input size. (E.g. if you go to bars and stripes for 8x8)

**Second assignment b) MNIST:**

Do the 2nd assigment with MNIST. I recommend using a Google Colab file or any other platform with GPU support to train the MNIST. On a CPU it will take quite some time. For this you will have to set up a Colab file that can run all the codes.

Use the 'Load_MNIST.py' file from the 'RBM' folder to load the files and explore the data set. I put some basic functions to open and show the images. And also how we can blank certain parts of the image to later reconstruct it with an RBM. Get familiar with the dataset and learn how to load and handle it.

- MNIST is a bit more challenging, but also more interesting to play with.
- The file `Load_MNIST.py` shows you already how to load MNIST files and how to make them from grey scale to black and white. (RBMs only work for binary inputs) Try to understand this file and use it or something similar to load your data.
- We have to convert the 28x28 pixel images into a 784 dimensional vector and put them into a RBM with the same input dimension. Find out how to do this with for example 'numpy'.
- What happens if you train the RBM just on one type of numbers. For example extract all the number 3 images from the data set and just train on them. To do so, the `data/` folder contains also the labels of the training set to find all images that contain e.g. a 3.


**Third assignment a):**
- If you go to `RBM_helper.py` file and change function `draw_sample(self, sample_length)` such that the Gibbs sampling is not started with an random vector, but with an actual image from the training set. What happens then? (For this task you will have to change `RBM_helper.py`). Play with the number of Gibbs steps. Does it make a big difference?

**Third assignment b):**
- Use now the images from `Bars_and_Stripes.py` that are partially blanked, which are stored in the variable 'subset' and also saved as a numpy file `blanked_bars_and_stripes.npy` if you run `Bars_and_Stripes.py`. You can think of these images as partially damaged and we would like to reconstruct them with our trained RBM.
- For this you will have to make changes in the `RBM_helper.py` file. In the function `draw_sample(self, sample_length)` we so far used a random vector to start the Gibbs sampling and make `sample_length` Gibbs setps until we obtain an output. Now we would like to start the Gibbs sampling with the 'damaged' images and see if we can reconstruct them. Do that and show your results.

If the outcome of the reconstructed image is just a random sample from the dataset and not the image that was damaged, you will have to fix the pixels of the input image that are not damaged after every Gibbs step. Give it a try and check if your results improve.

**Fourth Assigment:**

Now that you understand the basics of an RBM we can start doing quantum physics with RBMs. Open `RBM_train_on_H2_data_EXERCISE.ipynb` and read the introduction there. In the folder `H2_data/` you will find sampled data from a $H_2$ molecule.

- You will have to understand how to calculate the energy from these samples.
- Train a RBM on the sampling data and sample new data from the trained RBM.
- Check the energy of the quantum state that is learned by the RBM and compare it to the original data.
- This [paper](https://arxiv.org/pdf/1512.06860.pdf) should explain th enecessary physics behind the data. Equation (1) shows the Hamiltonian (The energy function). Figure 3a) shows how the energy should behave with the radius. The coefficients in Equation (1) can be found in the appendix of the paper and also in the file `H2_data/H2_coefficients.txt`.
- We simplified the problem a bit by changing the signs of the coefficients $g_4$ and $g_5$. You will have to do this as well to obtain the correct energy.


# Google Colab File
For those who cannot make all the packages run on their computer, they can create a google Colab file. For this they need a google / gmail account and access their google drive.

The Homework for fully connected NN and CNN should be done in the Colab file directly. The other homework can be done either on Colab or locally.


**Pytorch and Colab**
The following [link](https://drive.google.com/file/d/1Ppzwrb8yCDe7yF_OXJql-ygCeVOKxcs4/view?usp=sharing) provides a Colab file that helps you set up everything on Colab for pytorch and shows you how to load files on Colab.

Don't forget to activate the GPUs/TPUs (Go to Edit and then notebook settings) when using keras or pytorch. Beware: TPUs work only with tensorflow and keras.
