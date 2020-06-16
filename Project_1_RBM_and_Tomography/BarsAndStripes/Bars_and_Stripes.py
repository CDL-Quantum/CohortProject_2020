import numpy as np
import matplotlib.pyplot as plt

dataset = []
set_size = 10000

for k in range(set_size):
    init = np.zeros((4,4))
    nr_of_stripes = np.random.randint(low = 1, high= 4, size = 1)
    columns = np.random.randint(4, size = nr_of_stripes)

    for i in columns:
        init[:, i] = 1.0

    dataset.append(init)
    dataset.append(init.T)   
    
np.save('bars_and_stripes', dataset)   

subset = dataset[0:10]
for img in subset:
    img[0:2,0:2]=0

np.save('blanked_bars_and_stripes', subset) 
# plt.imshow(dataset[0])
# plt.show()  

# plt.imshow(subset[0])
# plt.show()