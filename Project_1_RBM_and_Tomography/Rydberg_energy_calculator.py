import numpy as np
import torch

def Vij(R, N):

    vij = torch.zeros(N, N)
    
    for i in range(N-1):
        for j in range(i+1,N):
        
            vij[i][j] = 1./(R*(j-i))**6.

    return vij 

    #vij = torch.zeros(N-1)
    #for k in range(0, N-1):
    #    vij[k] = 1./(R*(k+1))**6.

    #return vij
     
def target_psi(samples):
    # only for N = 2
    coeffs = torch.zeros(samples.shape[0])
    for i in range(samples.shape[0]):
        if torch.equal(samples[i,:], torch.tensor([0,0])):
            coeffs[i] = 0.973349
        if torch.equal(samples[i,:], torch.tensor([1,0])):
            coeffs[i] = 0.159871
        if torch.equal(samples[i,:], torch.tensor([0,1])):
            coeffs[i] = 0.159871
        if torch.equal(samples[i,:], torch.tensor([1,1])):
            coeffs[i] = 0.0383912
    return coeffs
 
def energy(samples, RBM_wvfn):

    # hardset Rydberg Hamiltonian parameters
    R = 1
    h = 1
    Omega = 1
    
    samples = samples.type(torch.int64)

    N = samples.shape[1]
    vij = Vij(R, N)

    # to +/- 1
    samples = samples*2 - 1

    diagonal1 = 0
    diagonal2 = 0
    # pair-wise term
    for i in range(N-1):
        for j in range(i+1,N):
            
            diagonal1 = diagonal1 - vij[i][j] * (
                samples[:,i] * samples[:,j] + samples[:,i] + samples[:,j]
            )
            #diagonal1 = diagonal1 - vij[j-i-1] * (
            #    samples[:,i] * samples[:,j] + samples[:,i] + samples[:,j]
            #)
            
        # longitudinal field term
        diagonal2 = diagonal2 - Omega * samples[:,i]

    # don't get the Nth term in the above loop
    diagonal2 = diagonal2 - Omega * samples[:,N-1]

    # back to 0,1
    samples = (samples + 1) / 2

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

    return (diagonal1 + diagonal2 + off_diagonal).sum() / (N*samples.shape[0])
