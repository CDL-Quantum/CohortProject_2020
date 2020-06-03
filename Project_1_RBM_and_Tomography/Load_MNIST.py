# Example file how to load MNIST from local folder.
# Download files from this site http://yann.lecun.com/exdb/mnist/
# Or use the ones in the folder RBM/data


import struct
import gzip
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import urllib.request
import os

def read_idx(filename):
    with gzip.open(filename, 'rb') as f:
        zero, data_type, dims = struct.unpack('>HBB', f.read(4))
        shape = tuple(struct.unpack('>I', f.read(4))[0] for d in range(dims))
        return np.fromstring(f.read(), dtype=np.uint8).reshape(shape)


directory = "data/"
if not os.path.exists(directory):
    os.makedirs(directory)

filename = "train-images-idx3-ubyte.gz"
if not os.path.exists(directory+filename):
    url = 'http://yann.lecun.com/exdb/mnist/train-images-idx3-ubyte.gz'
    urllib.request.urlretrieve(url, directory+filename)

# path to the .gz file
path = directory + filename
# Run the loader function
images = read_idx(path)
images = images/255 # normalize grey scale values between 0 and 255 to values between 0 and 1
images = np.where(images > 0.5, 1, 0) # This line makes the images black and white instead of grey scale. This is important, because the RBM is only binary and we cannot feed values different than 0 and 1.

# pick first 10 images of the whole dataset and make the upper
# right corner invisible for the reconstruction
subset = images[0:10]
for img in subset:
    img[0:15,0:15]=0

# Show the images of this subset
fig=plt.figure(figsize=(2, 5))
columns = 2
rows = 5
for i in range(0, columns*rows):
    fig.add_subplot(rows, columns, i+1)
    plt.imshow(subset[i])
plt.show()

# This is just another function how we can look at them without the axis
f, axarr = plt.subplots(5, 2)
kk = 0
for i in range(0,5):
    for j in range(0,2):
        axarr[i,j].imshow(subset[kk], cmap = cm.Greys_r)
        # axarr[i,j].set_title('title')
        axarr[i,j].axis('off')
        kk +=1
plt.show()

# Save these images to use them later as a test set for your RBM.
np.save(directory + 'RBM_blanked_test_images', subset)
