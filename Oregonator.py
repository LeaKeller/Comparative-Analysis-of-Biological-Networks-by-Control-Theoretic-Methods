# ---------------------------------------------------------------------------------------------------
# Trials made with the original reaction system implemented with BioNetGen that is not consequent. 
# The trials are left in comments. (in the other code, to complete here)
# ---------------------------------------------------------------------------------------------------

#import bionetgen
#from pysb.importers import bngl
#from pysb.core import ComponentSet
#from pysb.integrate import odesolve
#from pysb.pathfinder import set_path
#set_path("bng", "./code-path-folder/rl-env/Lib/site-packages/bionetgen/bng-win")

# Model (without bindings sites)
# Load the Oregonator.bngl env and assign it to the env variable
#env = bngl.model_from_bngl("Oregonator.bngl")
#run_env = bionetgen.run("Oregonator.bngl")

# Trial by reading the data from the BioNetGen (BNGL) .cdat file (not consequent)
#data = np.loadtxt("Oregonator.cdat", comments="#")

# ------------------------------------------------------------------------------------------------------------------------------------
# To pursue : 
# Implement the policy and algorithm
# Try the version where f is the optimal factor to be found among different combinations of concentrations/kinetics possibilities
# ------------------------------------------------------------------------------------------------------------------------------------

from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

# Implemented methods
methods = ['TDLearning', 'MonteCarlo']

# Define the Oregonator environment
class OregonatorEnv:      
    def __init__(self, oregonator): 
        """ 
        Constructor of the environment Oregonator.
        """
        super().__init__()

        self.oregonator = oregonator
        self.states = self.__states()
        print(self.states)
        for i in range(0,2):
            self.n_states = len(self.states[i])
        self.actions = self.__actions()
        #self.n_actions = len(self.actions)
        #self.rewards = self.__rewards() # in the actions for now

    def __oregonator(states, t):
        k1 = 1.34  # M^-1 sec^-1
        k2 = 1.6e9  # M^-1 sec^-1
        k3 = 8e3  # M^-1 sec^-1
        k4 = 4e7  # M^-1 sec^-1
        k5 = 1  # sec^-1
        n = 20
        f = np.random.randint(0, n) # CAREFUL since it is changing at integration level! Needs to be a fix value each time!
        
        A = 0.06
        B = 0.02
        
        q = 8.375e-6  
        
        X, Y, Z = states
    
        epsilon1 = k5/(k3*A)
        epsilon2 = (2*k4*k5)/(k2*k3*A)
        omega = 1
    
        dXdt = 1/epsilon1*(q*Y-X*Y+(X-X**2))
        dYdt = 1/epsilon2*(-q*Z-X*Y+f*Z)
        dZdt = omega*(X-Z)
    
        return dXdt, dYdt, dZdt
    
    def __states(self, states, t, f_action):
        self.tfinal = 300
        self.state_01 = 0.5
        self.state_02 = 0.1
        self.state_03 = 0.2
        t_eval = np.linspace(0, self.tfinal, 2000)

        X, Y, Z = states
        
        # Initial conditions
        self.state_0 = [self.state_01, self.state_02, self.state_03]
        t_eval = np.linspace(0, 300, 2000)
        
        # Solve ODE
        sol = odeint(self.__oregonator, self.state_0, t_eval)
        X, Y, Z = sol.T
        
        return X, Y, Z    

    def __next_states(self):        
        state = self.__states
        next_states = state[:,-1]
        return next_states
    
    def __unstable_regime_states(self, states, t, f_action):
        self.tfinal = 300
        self.state_01 = 0.5
        self.state_02 = 0.1
        self.state_03 = 0.2
        t_eval = np.linspace(0, self.tfinal, 2000)

        X, Y, Z = states
        
        # Initial conditions
        self.state_0 = [self.state_01, self.state_02, self.state_03]
        t_eval = np.linspace(0, 300, 2000)
        
        # Solve ODE
        sol = odeint(self.__unstable_regime, self.state_0, t_eval)
        X, Y, Z = sol.T
        
        return X, Y, Z  

    def __actions(self):
        def __unstable_regime(self, states, t):
            # Define constants as per the model
            k1 = 1.34  # M^-1 sec^-1
            k2 = 1.6e9  # M^-1 sec^-1
            k3 = 8e3  # M^-1 sec^-1
            k4 = 4e7  # M^-1 sec^-1
            k5 = 1  # sec^-1
            f_unstable = np.random.randint(1, 1+np.sqrt(2))

            A = 0.06 
            B = 0.02 
            
            q = 8.375e-6  
            
            X, Y, Z = states
        
            epsilon1 = k5/(k3*A)
            epsilon2 = (2*k4*k5)/(k2*k3*A)
            omega = 1
        
            dXdt = 1/epsilon1*(q*Y-X*Y+(X-X**2))
            dYdt = 1/epsilon2*(-q*Z-X*Y+f_unstable*Z)
            dZdt = omega*(X-Z)
        
            return dXdt, dYdt, dZdt
    
        self.tfinal = 300
        self.state_01 = 0.5
        self.state_02 = 0.1
        self.state_03 = 0.2
        t_eval = np.linspace(0, self.tfinal, 2000)
        
        # Initial conditions
        self.state_0 = [self.state_01, self.state_02, self.state_03]
        t_eval = np.linspace(0, 300, 2000)
        
        # Solve ODE
        sol = odeint(self.__unstable_regime, self.state_0, t_eval)
        X, Y, Z = sol.T
        
        self.state_target = self.__unstable_regime()
        ##Take and action if the system is stable
        #if state_target == ...
        #   reward = self.__states  # continue (the reward is set up as the state for now)
        #else:
        #    print("The system is stable")
        #return reward  

# --------------------------------------------------
# Simple plots of the synamics and phase without RL
# --------------------------------------------------

def oregonator_function(x, t):
    # Define constants as per the model
    k1 = 1.34  # M^-1 sec^-1
    k2 = 1.6e9  # M^-1 sec^-1
    k3 = 8e3  # M^-1 sec^-1
    k4 = 4e7  # M^-1 sec^-1
    k5 = 1  # sec^-1
    n = 2
    f = 0.5 # CAREFUL since it is changing at integration level! Needs to be a fix value!
    
    A = 0.06
    B = 0.02
    
    q = 8.375e-6  
    
    X, Y, Z = x

    epsilon1 = k5/(k3*A)
    epsilon2 = (2*k4*k5)/(k2*k3*A)
    omega = 1

    dXdt = 1/epsilon1*(q*Y-X*Y+(X-X**2))
    dYdt = 1/epsilon2*(-q*Z-X*Y+f*Z)
    dZdt = omega*(X-Z)
    
    return dXdt, dYdt, dZdt

def draw_oregonator_dynamics():
    tfinal = 300
    x01 = 0.5
    x02 = 0.1
    x03 = 0.2
    t_span = (0, 300)
    t_eval = np.linspace(0, tfinal, 2000)
    
    # Initial conditions
    x0 = [x01, x02, x03]
    t_eval = np.linspace(0, 300, 2000)
    
    # Solve ODE
    sol = odeint(oregonator_function, x0, t_eval)
    
    X, Y, Z = sol.T

    t = t_eval
    
    tau = t * 6.21  # Time scaling factor
    
    plt.figure(figsize=(10,6))
    plt.plot(tau, np.log10(X) + 10.30)
    plt.grid(True)
    plt.title("log[HBrO2] vs Time")
    plt.xlabel("Time (sec)")
    plt.ylabel("log[HBrO2]")
    
    plt.figure(figsize=(10,6))
    plt.plot(tau, np.log10(Y) + 6.52)
    plt.grid(True)
    plt.title("log[Br^-] vs Time")
    plt.xlabel("Time (sec)")
    plt.ylabel("log[Br^-]")
    
    plt.figure(figsize=(10,6))
    plt.plot(tau, np.log10(Z) + 7.62)
    plt.grid(True)
    plt.title("log[Ce(IV)] vs Time")
    plt.xlabel("Time (sec)")
    plt.ylabel("log[Ce(IV)]")

    plt.figure()
    plt.plot(np.log10(Y) + 6.5228, np.log10(X) + 10.2988)
    plt.grid(True)
    plt.title("Phase plot of log[HBrO2] vs log[Br^-] for f=0.5")
    plt.xlabel("log[Br^-]")
    plt.ylabel("log[HBrO2]")

    plt.figure()
    plt.plot(np.log10(Y) + 6.5228, np.log10(Z) + 7.62)
    plt.grid(True)
    plt.title("Phase plot of log[Ce(IV)] vs log[Br^-] for f=0.5")
    plt.xlabel("log[Br^-]")
    plt.ylabel("log[Ce(IV)]")
    
    plt.show()


if __name__ == "__main__":
    draw_oregonator_dynamics()

    # ------------------
    # Environment setup
    # ------------------
    def oregonator_function(x, t):
        k1 = 1.34  # M^-1 sec^-1
        k2 = 1.6e9  # M^-1 sec^-1
        k3 = 8e3  # M^-1 sec^-1
        k4 = 4e7  # M^-1 sec^-1
        k5 = 1  # sec^-1
        n = 20
        f = 1 #np.random.randint(0, n) # CAREFUL since it is changing at integration level! Needs to be a fix value!
        
        A = 0.06
        B = 0.02
        
        q = 8.375e-6  
        
        X, Y, Z = x
    
        epsilon1 = k5/(k3*A)
        epsilon2 = (2*k4*k5)/(k2*k3*A)
        omega = 1
    
        dXdt = 1/epsilon1*(q*Y-X*Y+(X-X**2))
        dYdt = 1/epsilon2*(-q*Z-X*Y+f*Z)
        dZdt = omega*(X-Z)
        
        return dXdt, dYdt, dZdt
    
    tfinal = 300
    x01 = 0.5
    x02 = 0.1
    x03 = 0.2
    t_span = (0, 300)
    t_eval = np.linspace(0, tfinal, 2000)
    
    # Initial conditions
    x0 = [x01, x02, x03]
    t_eval = np.linspace(0, 300, 2000)
    
    # Solve ODE
    sol = odeint(oregonator_function, x0, t_eval)
    
    X, Y, Z = sol.T

    t = t_eval
    
    tau = t * 6.21  # Time scaling factor

    oregonator = plt.figure(figsize=(10,6))
    ax = oregonator.add_subplot()
    ax.plot(np.log10(Y) + 6.5228, np.log10(X) + 10.2988)
    ax.grid(True)
    ax.set_title("Phase plot of log[HBrO2] vs log[Br^-] for f=1")
    ax.set_xlabel("log[Br^-]")
    ax.set_ylabel("log[HBrO2]")
    
    #env = OregonatorEnv(oregonator, states, actions)

    plt.show()