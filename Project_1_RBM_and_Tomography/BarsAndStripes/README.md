**Bars and Stripes:**
- Train the RBM on the Bars and Stripes data.
- In the file `Bars_and_Stripes.py` you can see how the data set is generated and saved to a `.npy` file.
- The bars and stripes "images" are only 4x4. Later you should scale them to bigger images.
- We have to convert the 4x4 pixel images into a 16 dimensional vector and put them into a RBM with the same input dimension. Find out how to do this with for example 'numpy'.
- Train the RBM (see the Task1 jupyter notebook for an example)
- Sample from the RBM and see what you obtain. Are the bars and stripes you sample good or are they blurry or just random noise? What happens if you increase the training epochs or the steps of the Gibbs sampling? Be aware that you have to transform the output of the RBM back to a 4x4 image. The same way you transformed the 4x4 image before to a 16 dimensional vector.
- Explore what happens if you change the number of hidden units. In the code so far we set the number of visible and hidden units equal. Can we still learn anything if we take e.g. 4 times less hidden units? How do the images look like? Explore and discuss! How does this scale with the input size. (E.g. if you go to bars and stripes for 8x8)