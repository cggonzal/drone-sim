import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import OdeSolution

def plot_results(sol: OdeSolution) -> None:
    # TODO plot results from sol
    return

def quad_model(t: np.ndarray, y: np.ndarray):
# reference https://sal.aalto.fi/publications/pdf-files/eluu11_public.pdf
    # TODO quad dynamics
    return   

time_start = 0
time_end = 10

# TODO define the state vector
#y0 = np.array([])

#sol = solve_ivp(quad_model, [time_start, time_end], y0)

# TODO
#plot_results(sol)
