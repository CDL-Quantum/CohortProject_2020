import numpy as np
import torch
from itertools import chain
from torch import nn
from torch.nn import functional as F
from torch.nn.utils import parameters_to_vector
from itertools import chain
from math import ceil

class RBM():

    def __init__(self, n_vis, n_hin):
        super(RBM, self).__init__()
        self.n_vis = n_vis
        self.n_hin = n_hin
        
        self.initialize_parameters()

    def initialize_parameters(self):
        self.weights = torch.randn(
            self.n_hin, 
            self.n_vis, 
            dtype=torch.double
        ) / np.sqrt(self.n_vis)
        
        self.visible_bias = torch.zeros(self.n_vis, dtype=torch.double)
        self.hidden_bias = torch.zeros(self.n_hin, dtype=torch.double)

    def effective_energy(self, v):
        v = v.to(self.weights)
        visible_bias_term = torch.matmul(v, self.visible_bias)
        hid_bias_term = F.softplus(F.linear(v, self.weights, self.hidden_bias)).sum(-1)

        return -(visible_bias_term + hid_bias_term)

    def effective_energy_gradient(self, v):
        v = (v.unsqueeze(0) if v.dim() < 2 else v).to(self.weights)
        prob = self.prob_h_given_v(v)

        W_grad = -torch.matmul(prob.transpose(0, -1), v)
        vb_grad = -torch.sum(v, 0)
        hb_grad = -torch.sum(prob, 0)
        return W_grad, vb_grad, hb_grad

    def prob_v_given_h(self, h):
        return (
            torch.matmul(h, self.weights.data, out=None)
            .add_(self.visible_bias.data)
            .sigmoid_()
            .clamp_(min=0, max=1)
        )

    def prob_h_given_v(self, v):
        return (
            torch.matmul(v, self.weights.data.t(), out=None)
            .add_(self.hidden_bias.data)
            .sigmoid_()
            .clamp_(min=0, max=1)
        )

    def sample_v_given_h(self, h):
        v = self.prob_v_given_h(h)
        v = torch.bernoulli(v)  # overwrite v with its sample
        return v

    def sample_h_given_v(self, v):
        h = self.prob_h_given_v(v)
        h = torch.bernoulli(h)  # overwrite h with its sample
        return h

    def draw_samples(self, k, initial_state):
        v = (initial_state.clone()).to(self.weights)
        h = torch.zeros(*v.shape[:-1], self.n_hin).to(self.weights)

        for _ in range(k):
            h = self.sample_h_given_v(v)
            v = self.sample_v_given_h(h)

        return v

    def wavefunction(self, v):
        return (-self.effective_energy(v)).exp().sqrt()
    
    def gradients(self, batch):
        grads_W, grads_vb, grads_hb = self.effective_energy_gradient(batch)
        grads_W /= float(batch.shape[0])
        grads_vb /= float(batch.shape[0])
        grads_hb /= float(batch.shape[0])

        return grads_W, grads_vb, grads_hb

    def compute_batch_gradients(self, k, pos_phase_batch, neg_phase_batch):
        grads_W, grads_vb, grads_hb = self.gradients(pos_phase_batch)

        vk = self.draw_samples(k, neg_phase_batch)
        neg_grads_W, neg_grads_vb, neg_grads_hb = self.gradients(vk)

        grads_W -= neg_grads_W
        grads_vb -= neg_grads_vb
        grads_hb -= neg_grads_hb

        return grads_W, grads_vb, grads_hb
  
    def shuffle_data(self, data, batch_size):
        permutation = torch.randperm(data.shape[0])
        data = [
            data[batch_start : (batch_start + batch_size)]
            for batch_start in range(0, len(data), batch_size)
        ]
        return data

    def params(self):
        return self.weights, self.visible_bias, self.hidden_bias

    def update_params(self, grads, lr):
        self.weights -= lr*grads[0]
        self.visible_bias -= lr*grads[1]
        self.hidden_bias -= lr*grads[2]

    def train(self, input_data, k=10, batch_size=100, lr=0.01):
          
        num_batches = ceil(input_data.shape[0] / batch_size)
        pos_batches = self.shuffle_data(input_data, batch_size)
        neg_batches = self.shuffle_data(input_data, batch_size)

        for b in range(num_batches):
            all_gradients = self.compute_batch_gradients(k, pos_batches[b], neg_batches[b])
            
            self.update_params(all_gradients, lr)

    def partition_function(self, space):
        logZ = (-self.effective_energy(space)).logsumexp(0)
        return logZ.exp()

    def generate_hilbert_space(self):
        dim = np.arange(2 ** self.n_vis)
        space = ((dim[:, None] & (1 << np.arange(self.n_vis))) > 0)[:, ::-1]
        space = space.astype(int)
        return torch.tensor(space, dtype=torch.double)

    def psi(self):
        space = self.generate_hilbert_space()
        return self.wavefunction(space) / self.partition_function(space).sqrt()
