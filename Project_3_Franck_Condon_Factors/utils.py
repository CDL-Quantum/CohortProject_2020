import warnings
import numpy as np
import strawberryfields as sf
import matplotlib.pyplot as plt

def get_state(
    t: np.ndarray,
    U1: np.ndarray,
    r: np.ndarray,
    U2: np.ndarray,
    alpha: np.ndarray,
    backend: str = "gaussian",
    backend_options: dict = {},
    loss: float = 0.0,
    ) -> list:

    if not 0 <= loss <= 1:
        raise ValueError("Loss parameter must take a value between zero and one")

    n_modes = len(t)

    eng = sf.Engine(backend=backend, backend_options=backend_options)

    if np.any(t != 0):
        prog = sf.Program(n_modes * 2)
    else:
        prog = sf.Program(n_modes)

    # pylint: disable=expression-not-assigned,pointless-statement
    with prog.context as q:

        if np.any(t != 0):
            for i in range(n_modes):
                sf.ops.S2gate(t[i]) | (q[i], q[i + n_modes])

        sf.ops.Interferometer(U1) | q[:n_modes]

        for i in range(n_modes):
            sf.ops.Sgate(r[i]) | q[i]

        sf.ops.Interferometer(U2) | q[:n_modes]

        for i in range(n_modes):
            sf.ops.Dgate(np.abs(alpha[i]), np.angle(alpha[i])) | q[i]

        if loss:
            for _q in q:
                sf.ops.LossChannel(1 - loss) | _q

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning, message="Cannot simulate non-")

    return eng.run(prog).state


def lorentzian_smoothing(sticks, gamma=100.0, points=1000):
    """
    Lorentzian smothing of sticks containting intensities of energies. Intensity here is a 
    broad term, it can be the FCF factors, counts, etc.
    
    Args:
        sticks: 2D numpy array [[energies], [intensities]]
        points: Number of points to divide the interval into
        gamma: The width (at half maximum) of Lorentzian function
        
    Return:
        A 2D numpy array of size points*points containing [[energies], [smoothed_intensities]]
    """
    energies = sticks[0, :]
    intensities = sticks[1, :]
    
    _min, _max = min(energies), max(energies)
    emin = _min - 0.1 * (_max - _min)
    emax = _max + 0.1 * (_max - _min) 
    
    x = np.linspace(emin, emax, points)
    y = np.zeros(points)
    for energy, intensity in zip(energies, intensities):
        y += intensity / (1 + ((x - energy) **2 / ((gamma / 2) **2)))
        
    return np.vstack([x, y])


def plot_spectrum(sticks, description, xlim=None, smoothing=True, gamma=100.0, points=None):
    """
    Args:
        sticks: Spectrum sticks to plot. 2D numpy array [[energies], [intensities]]
        xlim: tuple (xmin, xmax)
        description: dictionary, containing title, xlabel, ylabel
        smoothing: If True, does Lorentzian smoothing and plots it as well
        points: If smoothing=True, indicates the number of points to use within the energy spaning interval
    """
    energies = sticks[0, :]
    FCFs = sticks[1, :]
    
    if xlim:
        plt.xlim(xlim[0], xlim[1])
        
    plt.vlines(energies, 0, FCFs)
    if (smoothing):
        if not points:
            points = int(max(energies) - min(energies))
        ls = lorentzian_smoothing(sticks, gamma=gamma, points=points)
        plt.plot(ls[0, :], ls[1, :])
    
    plt.title(description.get('title', ''))
    plt.xlabel(description.get('xlabel', ''))
    plt.ylabel(description.get('ylabel', ''))
    plt.show()
    

def plot_spectrum_from_samples(energies, description, xlim=None, smoothing=True, gamma=100.0,  points=1000): 
    energy_counts = {}
    for e in energies:
        energy_counts[e] = energy_counts.get(e, 0) + 1
    
    sticks = np.vstack([list(energy_counts.keys()), list(energy_counts.values())])
    plot_spectrum(sticks, description, xlim, smoothing, gamma, points)
