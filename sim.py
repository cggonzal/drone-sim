import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import OdeSolution

def plot_results(sol: OdeSolution) -> None:
    # TODO plot results from sol
    return

def get_epsilon_dot_dot() -> np.ndarray:
    pass

def get_n_dot_dot(X: np.ndarray, u: np.ndarray) -> np.ndarray:
# X = [x, y, z, x_dot, y_dot, z_dot, phi, theta, psi, phi_dot, theta_dot, psi_dot]   
    g = -9.81
    m = 1 # NOTE change mass as needed
    T = u[0 ,0]

    phi = X[6]
    theta = X[7]
    psi = X[8]

    z_dot_dot = g + (T / m) * np.cos(theta) * np.cos(phi)
    
    n_dot_dot = np.array([x_dot_dot, y_dot_dot, z_dot_dot])
    return n_dot_dot

def quad_model(t: np.ndarray, X: np.ndarray): # TODO add control input "u" as parameter to quad_model
# reference https://sal.aalto.fi/publications/pdf-files/eluu11_public.pdf
# input state: 
# X = [x, y, z, x_dot, y_dot, z_dot, phi, theta, psi, phi_dot, theta_dot, psi_dot]

#output state: 
#X_dot = [x_dot, y_dot, z_dot, x_dot_dot, y_dot_dot, z_dot_dot, phi_dot, theta_dot, psi_dot, phi_dot_dot, theta_dot_dot, psi_dot_dot]
    x = X[0]
    y = X[1]
    z = X[2]
    x_dot = X[3]
    y_dot = X[4]
    z_dot = X[5]
    phi = X[6]
    theta = X[7]
    psi = X[8]
    phi_dot = X[9]
    theta_dot = X[10]
    psi_dot = X[11]
    
    epsilon_dot_dot = get_epsilon_dot_dot()
    n_dot_dot = get_n_dot_dot()
    
    X_dot = np.zeros((12, 1))
    X_dot[0] = x_dot
    X_dot[1] = y_dot
    X_dot[2] = z_dot
    X_dot[3] = epsilon_dot_dot[0, 0]
    X_dot[4] = epsilon_dot_dot[1, 0]
    X_dot[5] = epsilon_dot_dot[2, 0]
    X_dot[6] = phi_dot
    X_dot[7] = theta_dot
    X_dot[8] = psi_dot
    X_dot[9] = n_dot_dot[0, 0]
    X_dot[10] = n_dot_dot[1, 0]
    X_dot[11] = n_dot_dot[2, 0]

    return X_dot 

time_start = 0
time_end = 10

# TODO define the initial state vector
y0 = np.array([])

sol = solve_ivp(quad_model, [time_start, time_end], y0)

# TODO
#plot_results(sol)
