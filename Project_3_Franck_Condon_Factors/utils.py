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


def lorentzian_smoothing(sticks, gamma=100.0):
    energies = sticks[0, :]
    FCFs = sticks[1, :]
    
    emin = min(energies) - 1000
    emax = max(energies) + 1000
    
    x = np.linspace(emin, emax, int(emax - emin))
    y = 0
    for energy, FCF in zip(energies, FCFs):
        y += FCF / (1 + ((x - energy) **2 / ((gamma / 2) **2)))
        
    return np.vstack([x, y])


def plot_spectrum(sticks, title, smoothing=True):
    energies = sticks[0, :]
    FCFs = sticks[1, :]
    
    plt.vlines(energies, 0, FCFs)
    if (smoothing):
        ls = lorentzian_smoothing(sticks)
        plt.plot(ls[0, :], ls[1, :])
    
    plt.title(title)
    plt.xlabel('Energy ($cm^{-1}$)')
    plt.ylabel('Intensity')
    plt.show()
    

def plot_spectrum_from_samples(energies, title, smoothing=True): 
    energy_counts = {}
    for e in energies:
        energy_counts[e] = energy_counts.get(e, 0) + 1
    
    sticks = np.vstack([list(energy_counts.keys()), list(energy_counts.values())])
    plot_spectrum(sticks, title, smoothing)
