**MNIST:**

I recommend using a Google Colab file or any other platform with GPU support to train on MNIST data. On a CPU it will take quite some time. For this you will have to set up a Colab file that can run all the codes.

Use the 'Load_MNIST.py' file to load the files and explore the data set. I put some basic functions to open and show the images. And also how we can blank certain parts of the image to later reconstruct it with an RBM. Get familiar with the dataset and learn how to load and handle it.

- MNIST is a bit more challenging, but also more interesting to play with.
- The file `Load_MNIST.py` shows you already how to load MNIST files and how to make them from grey scale to black and white. (RBMs only work for binary inputs) Try to understand this file and use it or something similar to load your data.
- We have to convert the 28x28 pixel images into a 784 dimensional vector and put them into a RBM with the same input dimension. Find out how to do this with for example 'numpy'.
- What happens if you train the RBM just on one type of numbers. For example extract all the number 3 images from the data set and just train on them. To do so, the `data/` folder contains also the labels of the training set to find all images that contain e.g. a 3.

**Google Colab File** 

For those who cannot make all the packages run on their computer, they can create a google Colab file. For this they need a google / gmail account and access their google drive.

The Homework for fully connected NN and CNN should be done in the Colab file directly. The other homework can be done either on Colab or locally.

**Pytorch and Colab**

The following [link](https://drive.google.com/file/d/1Ppzwrb8yCDe7yF_OXJql-ygCeVOKxcs4/view?usp=sharing) provides a Colab file that helps you set up everything on Colab for pytorch and the how to load files on Colab.

Don't forget to activate the GPUs/TPUs (Go to Edit and then notebook settings) when using keras or pytorch. Beware: TPUs work only with tensorflow and keras.