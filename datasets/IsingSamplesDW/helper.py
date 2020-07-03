import numpy as np


def load_dataset(filename):
    """
    A function that loads the data set of samples of an Ising Hamiltonian given the filename.

    Args:
        filename (str): name of the file

    Returns:
        train (numpy.ndarray dtype=np.uint8): The training set with shape=(40000, 2041)
        val (numpy.ndarray dtype=np.uint8): The validation set with shape=(10000, 2041)
    """
    try:
        dataset = np.load(filename)
        train = dataset['train']
        val = dataset['validation']
        assert(train.shape == (40000, 2041))
        assert(val.shape == (10000, 2041))
        assert(all(np.unique(train) == [0, 1]))
        assert(all(np.unique(val) == [0, 1]))
        return train, val
    except IOError:
        raise IOError("make sure you point to the right dataset_{hashcode}.npz file")
    except AssertionError:
        raise AssertionError("The properties of the file are not correct")
    except KeyError:
        raise KeyError("Make sure the fields train and validation exist in the dataset")
