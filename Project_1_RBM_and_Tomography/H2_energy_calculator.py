import numpy as np
import torch

def find_sample_frequency(samples):
    freq = []
    basis = np.array([[0.0, 0.0], [0.0, 1.0], [1.0, 0.0], [1.0, 1.0]])
    basis = torch.from_numpy(basis)   
 
    for b in basis:
        l = 0
        for sample in samples:
            if (sample == b).all():
                l += 1
        freq.append(l)
    return np.array(freq) / len(samples)

def energy_from_freq(samples, coeff):
    f = find_sample_frequency(samples)
    c1, c2, c3, c4, c5, c6 = coeff[1:-1]
    
    diag = c1*one(f) + c2*Z0(f) + c3*Z1(f) + c4*Z0Z1(f) 
    off_diag = -c5*(2*np.sqrt(f)[[1,2]].prod() + 2*np.sqrt(f)[[0,3]].prod()) - c5*(2*np.sqrt(f)[[1,2]].prod() - 2*np.sqrt(f)[[0,3]].prod())

    return diag + off_diag

def one(freq):
    return freq.sum()

def Z0(freq):
    return freq[:2].sum() - freq[2:4].sum()

def Z1(freq):
    return freq[[0,2]].sum() - freq[[1,3]].sum()

def Z0Z1(freq):
    return freq[[0,3]].sum() - freq[[1,2]].sum()

def X0X1(wvfn, samples):
    samples = samples.type(torch.double)
    denom = wvfn(samples)
    samples = samples.type(torch.int64)
    tmp = (samples ^ 1).type(torch.double)
    num = wvfn(tmp)
    ratio = num / denom

    # matrix elements are +1
    return ratio.sum() / samples.size()[0]

def Y0Y1(wvfn, samples):
    denom = wvfn(samples)
    samples = samples.type(torch.int64)
    tmp = (samples ^ 1).type(torch.double)
    num = wvfn(tmp)
    ratio = num / denom
   
    # matrix elements are +/-1
    # -1: if sum(Sz) = 0 (mod 2)
    # +1: otherwise 
    parity = (samples.sum(dim=1)) % 2
    ratio *= (-1)**(parity+1)
    freq = find_sample_frequency(samples)
    return ratio.sum() / samples.size()[0] 

def energy(samples, coeff, RBM_wvfn):
    '''
    The Hamiltonian for H2, under a Bravyi-Kitaev transformation, is

    H = c1 + c2*\sigma_1^z + c3*\sigma_2^z + c4*\sigma_1^z \sigma_2^z +
        c5*\sigma_1^x \sigma_2^x + c5*\sigma_1^y \sigma_2^y

    (equation 37 in Xia et al)

    However, there's a sign problem in this Hamiltonian which complicates the
    learning problem quite a bit (the ground state wavefunction has sign 
    structure, so we require a different RBM and data in other bases!). 

    The sign problem can be alleviated by transforming the Hamiltonian with the
    unitary kron(1, \sigma^z). This manifests in flipping the sign of c5.
    '''
    

    # only need freqs for first four terms
    f = find_sample_frequency(samples)
    c1, c2, c3, c4, c5, c6 = coeff[1:-1]
    diagonal_part = c1*one(f) + c2*Z0(f) + c3*Z1(f) + c4*Z0Z1(f)

    # need wavefunction coefficients (unnorm) for other terms (off diag)
    offdiagonal_part = -c5*X0X1(RBM_wvfn, samples) - c5*Y0Y1(RBM_wvfn, samples)

    return diagonal_part + offdiagonal_part
