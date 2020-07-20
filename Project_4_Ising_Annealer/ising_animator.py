from matplotlib import animation, rc
import matplotlib.pyplot as plt
from IPython.display import HTML

class IsingAnimator():
    def __init__(self, system):
        self.system = system
        self.image = None
        
        # set up the figure
        self.fig, self.ax = plt.subplots()
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        
    def get_spins(self):
        return self.system.spins.reshape(-1, self.system.spins.shape[-1]).copy()

    def initialize_system(self):
        """Initializes the animation for use with matplotlib.animation.FuncAnimation()"""
        self.image = self.ax.imshow(self.get_spins(), 
                                    interpolation='none', cmap='gray')
        return (self.image,)

    def animate_system(self, T):
        """Performs a frame update of the current animation. 

        `T` is the temperature at which to perform the underlying MC update
        """
        E = self.system.mc_step(T)
        self.image.set_data(self.get_spins())

        self.ax.set_title("$T = {}$\n$E = {}$".format(T, E))

        return (self.image,)


    def run_animation(self, T_sched):
        """Creates an animation of the Ising Lattice simulation following the given annealing schedule

        `T_sched` is an annealing schedule, should be a list/array
        """

        anim = animation.FuncAnimation(self.fig, self.animate_system, init_func=self.initialize_system,
                                       frames=T_sched, interval=100, 
                                       blit=False, cache_frame_data=False)
        return HTML(anim.to_html5_video())