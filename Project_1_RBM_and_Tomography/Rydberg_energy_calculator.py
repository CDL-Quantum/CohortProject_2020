import numpy as np
import torch

def Vij(R,i,j):
    return 1./(R*abs(i-j))**6.

def energy(samples, RBM_wvfn):

    # hardset Rydberg Hamiltonian parameters
    R = 1.
    h = 1.
    Omega = 1.
    
    samples = samples.type(torch.int64)
    N = samples.shape[1]

    # to +/- 1
    samples = -(samples*2 - 1)

    diagonal = 0
    # pair-wise term
    for i in range(N-1):
        diagonal -= Omega*samples[:,i]
        for j in range(i+1,N):
            diagonal -= Vij(R,i,j)*(samples[:,i]*samples[:,j] + samples[:,i] + samples[:,j])
            
    # don't get the Nth term in the above loop
    diagonal -= Omega*samples[:,N-1]

    # back to 0,1
    samples = (-samples + 1) / 2

    # transverse-field term
    off_diagonal = 0
    for n in range(N):
        denom = RBM_wvfn(samples)
        # flip the nth qubit
        samples[:,n] = samples[:,n] ^ 1
        
        num = RBM_wvfn(samples)
        ratio = num / denom

        off_diagonal -= h*ratio
            
        # undo the flip
        samples[:,n] = samples[:,n] ^ 1

    return (diagonal + off_diagonal).sum() / (N*samples.shape[0])
