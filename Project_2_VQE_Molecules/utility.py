import math
import numpy as np
from openfermion import QubitOperator
from openfermion.hamiltonians import MolecularData
from openfermionpyscf import run_pyscf
from openfermion.transforms import get_fermion_operator, bravyi_kitaev, jordan_wigner

from tequila.grouping.binary_rep import BinaryHamiltonian
from tequila import QubitHamiltonian, quantumchemistry

def get_qubit_hamiltonian(mol, geometry, basis, charge=0, multiplicity=1, qubit_transf='bk'):
    '''
    Generating qubit hamiltonina of the given molecules with specified geometry, basis sets, charge and multiplicity
    Give its qubit form using specified transformation
    '''
    g = get_molecular_data(mol, geometry)
    mol = MolecularData(g, basis, multiplicity, charge)
    mol = run_pyscf(mol)

    ham = mol.get_molecular_hamiltonian()
    hamf = get_fermion_operator(ham)

    if qubit_transf == 'bk':
        hamq = bravyi_kitaev(hamf)
    elif qubit_transf == 'jw':
        hamq = jordan_wigner(hamf)
    else:
        raise(ValueError(qubit_transf, 'Unknown transformation specified'))

    return hamq

def convert_mol_data_to_xyz_format(mol_data):
    '''
    Convert nuclear geometry list to .xyz format.
    '''

    xyz_str = ''
    for atom in mol_data:
        xyz_str += atom[0] +' ' + ' '.join([f"{coord:.10f}" for coord in atom[1]]) +'\n'

    return xyz_str

def get_molecular_data(mol, geometry, xyz_format=False):
    '''
    Generate the molecular data of the specified molecule
    '''
    if mol == 'h2':
        mol_data = [
            ['H', [0, 0, 0]],
            ['H', [0, 0, geometry]]
        ]
    elif mol == 'lih':
        mol_data = [
            ['Li', [0, 0, 0]],
            ['H', [0, 0, geometry]]
        ]
    elif mol == 'h2o':
        # Giving symmetrically stretch H2O. ∠HOH = 107.6°
        # Geometry is distance between H-O
        angle = 107.6 / 2
        angle = math.radians(angle)
        x = geometry * math.sin(angle)
        y = geometry * math.cos(angle)
        mol_data = [
            ['O', [0, 0, 0]],
            ['H', [-x, y, 0]],
            ['H', [x, y, 0]]
        ]
    elif mol == 'n2':
        mol_data = [
            ['N', [0, 0, 0]],
            ['N', [0, 0, geometry]]
        ]
    elif mol == 'h4':
        mol_data = [
            ['H', [0, 0, 0]],
            ['H', [0, 0, geometry]],
            ['H', [0, geometry, 0]],
            ['H', [0, geometry, geometry]]
        ]
    elif mol == 'nh3':
        bondAngle = 107
        bondAngle = math.radians(bondAngle)
        cos = math.cos(bondAngle)
        sin = math.sin(bondAngle)

        # The idea is second and third vecctor dot product is cos(angle) * geometry^2.
        thirdyRatio = (cos - cos**2) / sin
        thirdxRatio = (1 - cos**2 - thirdyRatio**2) ** (1/2)
        mol_data = [
            ['H', [0.0, 0.0, geometry]],
            ['H', [0.0, sin * geometry, cos * geometry]],
            ['H', [thirdxRatio * geometry, thirdyRatio * geometry, cos * geometry]],
            ['N', [0.0, 0.0, 0.0]]
            ]

    else:
        raise(ValueError(mol, 'Unknown moleucles given'))

    if xyz_format:
        return convert_mol_data_to_xyz_format(mol_data)
    else:
        return mol_data


def get_number_qubit(H : QubitOperator):
    '''
    Return the number of qubits in H
    '''
    n_qub = -1
    for pw, val in H.terms.items():
        for ps in pw:
            n_qub = max(n_qub, ps[0])

    return n_qub + 1

def largest_first(commuting_graph_complement):
    '''
    Using a n x n binary matrix A where A[i, j] = 0 means the underlying ith item commutes with the
    jth item.
    Returns a dictionary whose values contains mutually commuting indices.
    '''
    n = commuting_graph_complement.shape[0]

    rows = commuting_graph_complement.sum(axis=0)
    ind = np.argsort(rows)[::-1]
    m = commuting_graph_complement[ind,:][:,ind]
    colors = dict()
    c = np.zeros(n, dtype=int)
    k = 0 #color

    for i in range(n):
        neighbors = np.argwhere(m[i,:])
        colors_available = set(np.arange(1, k+1)) - set(c[[x[0] for x in neighbors]])
        term = ind[i]
        if not colors_available:
            k += 1
            c[i] = k
            colors[c[i]] = [term]
        else:
            c[i] = min(list(colors_available))
            colors[c[i]].append(term)

    return colors

def pauli2binvec(pws, n):
    '''
    Turning list of pauli words into list of binary vectors in form [z, x]
    '''
    binvecs = []
    for pw in pws:
        vec = np.zeros(2*n)
        for ps in pw:
            qub = ps[0]
            w = ps[1]
            if w == 'Z':
                vec[qub] = 1
            elif w == 'X':
                vec[qub + n] = 1
            else:
                vec[qub] = 1
                vec[qub + n] = 1
        binvecs.append(vec)
    return binvecs

def anticommute(a, b):
    '''
    Return the binary symplectic inner product between two binary vectors a and b.

    Return: 0 or 1, commute or anti-commute.
    '''
    dim = len(a) // 2
    re = a[:dim] @ b[dim:] + b[:dim] @ a[dim:]

    return re % 2

def get_commuting_group(H : QubitOperator):
    '''
    Get a dictionary of mutually commuting groups with terms in Hp
    '''
    n = get_number_qubit(H)

    pws = []
    pws_val = []
    for pw, val in H.terms.items():
        pws.append(pw)
        pws_val.append(val)

    binvecs = pauli2binvec(pws, n)
    tnum = len(binvecs)
    comm_matrix = np.zeros((tnum, tnum))

    for i in range(tnum):
        for j in range(i+1, tnum):
            comm_matrix[i, j] = 1 - anticommute(binvecs[i], binvecs[j])

    comm_matrix = np.identity(tnum) + comm_matrix + comm_matrix.T
    colors = largest_first(1 - comm_matrix)

    dict = {}
    for key, indices in colors.items():
        dict[key] = QubitOperator.zero()
        for idx in indices:
            dict[key] += QubitOperator(term=pws[idx], coefficient=pws_val[idx])

    return dict

def get_qwc_unitary(H : QubitOperator):
    '''
    Get the unitary that transform commuting operators to qwc operators
    '''
    qh = QubitHamiltonian.from_openfermion(H)
    bh = BinaryHamiltonian.init_from_qubit_hamiltonian(qh)

    qwc, lag, sig = bh.single_qubit_form()

    num = len(lag)
    U = QubitOperator.identity()

    for idx in range(num):
        l = QubitHamiltonian.from_paulistrings(lag[idx].to_pauli_strings())
        s = QubitHamiltonian.from_paulistrings(sig[idx].to_pauli_strings())
        U *= 1 / 2 ** (1/2) * (l.to_openfermion() + s.to_openfermion())
    return U

def qubit_wise_commuting(a : QubitOperator, b : QubitOperator):
    '''
    Check if a and b are qubit-wise commuting.
    assume a and b have only one term
    '''
    ps_dict = {}

    pw, _ = a.terms.copy().popitem()

    for ps in pw:
        ps_dict[ps[0]] = ps[1]

    pw, _ = b.terms.copy().popitem()
    for ps in pw:
        if ps[0] in ps_dict:
            if ps[1] != ps_dict[ps[0]]:
                return False

    return True

def get_qwc_group(H : QubitOperator):
    '''
    Return a list of qubit-wise commuting fragments of H
    '''
    # Preparing all terms in H into a list
    qubit_ops = []
    for pw, val in H.terms.items():
        qubit_ops.append(QubitOperator(term=pw, coefficient=val))
    n = len(qubit_ops)

    # Making commutation matrix
    comm_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            comm_matrix[i, j] = qubit_wise_commuting(qubit_ops[i], qubit_ops[j])

    # Compute commuting fragments
    comm_matrix = np.identity(n) + comm_matrix + comm_matrix.T
    colors = largest_first(1 - comm_matrix)

    # Collect commuting fragments into a list of QubitOperators
    qwc_list = []
    qwc_list_idx = 0
    for key, indices in colors.items():
        qwc_list.append(QubitOperator.zero())
        for idx in indices:
            qwc_list[qwc_list_idx] += qubit_ops[idx]
        qwc_list_idx += 1
    return qwc_list

def obtain_PES(mol_data_list, basis, method):

    if method.lower() not in ['ccsd', 'cisd', 'fci', 'hf']:
        raise(ValueError("Method not recognized, implemented methods are 'ccsd', 'cisd', 'fci', 'hf'."))

    gridpoints = len(mol_data_list)

    energies = np.zeros(gridpoints)

    for i in range(gridpoints):

        moldata_R = quantumchemistry.Molecule(mol_data_list[i], basis)

        if method == 'cisd':
            result = moldata_R.compute_energy('detci', options={"detci__ex_level": 2})
        else:
            result = moldata_R.compute_energy(method)

        print("E = {} Eh".format(result))

        energies[i] = result

    return energies
