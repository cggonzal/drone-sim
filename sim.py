import numpy as np
from scipy.integrate import solve_ivp
from typing import Any
import matplotlib.pyplot as plt
#
def plot_results(sol: Any) -> None:
    final_t = sol.t[-1]
    num_points = 100

    t = np.linspace(0, final_t, num_points)
    z = sol.sol(t)
    
    plt.plot(t, z[2, :])    
    plt.xlabel("time(s)")
    plt.ylabel("z axis(m)")
    
    path = "/home/cgg/Desktop/figure.png"
    plt.savefig(path)
    print("saved figure to path:", path)

def get_W_minus_one_n_derivative(phi: float, theta: float, psi: float, 
        phi_dot: float, theta_dot: float, psi_dot: float) -> np.ndarray:
    # equation (12)
    zero_zero = 0
    zero_one = phi_dot * np.cos(phi) * np.tan(theta) + theta_dot * np.sin(phi) / (np.cos(theta)**2)
    zero_two = -phi_dot * np.sin(phi) * np.cos(theta) + theta_dot * np.cos(phi) / (np.cos(theta) ** 2)
    
    one_zero = 0
    one_one = -phi_dot * np.sin(phi)
    one_two = -phi_dot * np.cos(phi)
    
    two_zero = 0
    two_one = phi_dot * np.cos(phi) / np.cos(theta) + phi_dot * np.sin(phi) * np.tan(theta) / np.cos(theta)
    two_two = -phi_dot * np.sin(phi) / np.cos(theta) + theta_dot * np.cos(phi) * np.tan(theta) / np.cos(theta)
    
    W_minus_one_n_derivative = np.array([[zero_zero, zero_one, zero_two],
                                         [one_zero, one_one, one_two],
                                         [two_zero, two_one, two_two]])
    
    return W_minus_one_n_derivative

def get_W_n(phi:float, theta:float, psi: float) -> np.ndarray:
    # equation (4)
    W_n = np.array([[1, 0, -np.sin(theta)], 
                    [0, np.cos(phi), np.cos(theta) * np.sin(phi)],
                    [0, -np.sin(phi), np.cos(theta) * np.cos(phi)]])
    
    return W_n

def get_v(phi:float, theta:float, psi:float, phi_dot:float, theta_dot:float, psi_dot: float) -> np.ndarray:
    # equation (4)
    W_n = get_W_n(phi, theta, psi)
    
    n_dot = np.array([[phi_dot], [theta_dot], [psi_dot]])
    
    return W_n @ n_dot

def get_W_minus_one_n(phi:float, theta:float, psi:float) -> np.ndarray:
    # equation (4)
    W_minus_one_n = np.array([[1, np.sin(phi) * np.tan(theta), np.cos(phi) * np.tan(theta)],
                              [0, np.cos(phi), -np.sin(phi)],
                              [0, np.sin(phi) / np.cos(theta), np.cos(phi) / np.cos(theta)]])
    return W_minus_one_n

def get_u_from_omega(omega: np.ndarray, lift_constant: float, motor_torque_constant: float, 
                     arm_length: float) -> np.ndarray:
    # NOTE: counter-clockwise torque is positive, clockwise torque is negative
    # equations (7) and (8)
    motor_matrix_to_u_from_omega = np.array([[lift_constant, lift_constant, lift_constant, lift_constant],
                                             [0, -lift_constant * arm_length, 0, lift_constant * arm_length], 
                                             [-lift_constant * arm_length, 0, lift_constant * arm_length, 0], 
                                             [-motor_torque_constant, motor_torque_constant, -motor_torque_constant, motor_torque_constant]])

    return motor_matrix_to_u_from_omega @ (omega ** 2)

def get_omega_from_u(u: np.ndarray, lift_constant: float, motor_torque_constant: float, 
                     arm_length: float) -> np.ndarray:
    # inverts equations (7) and (8)
    motor_matrix_to_u_from_omega = np.array([[lift_constant, lift_constant, lift_constant, lift_constant],
                                             [0, -lift_constant * arm_length, 0, lift_constant * arm_length], 
                                             [-lift_constant * arm_length, 0, lift_constant * arm_length, 0], 
                                             [-motor_torque_constant, motor_torque_constant, -motor_torque_constant, motor_torque_constant]])
         
    motor_matrix_to_omega_from_u = np.linalg.inv(motor_matrix_to_u_from_omega)

    return np.sqrt(motor_matrix_to_omega_from_u @ u)


def get_v_dot(phi: float, theta: float, psi: float, phi_dot: float, theta_dot: float, 
        psi_dot: float, I: np.ndarray, u: np.ndarray, lift_constant: float, motor_torque_constant: float,
        arm_length: float) -> np.ndarray:
    # equation (11), assumes center of mass and body frame are aligned 
    I_xx = I[0, 0]
    I_yy = I[1, 1]
    I_zz = I[2, 2]

    tau_phi = u[1, 0]
    tau_theta = u[2, 0]
    tau_psi = u[3, 0]
    
    v = get_v(phi, theta, psi, phi_dot, theta_dot, psi_dot)
    p = v[0, 0]
    q = v[1, 0]
    r = v[2, 0]
    
    # NOTE gyroscopic force assumed to be zero 
    p_dot = (I_yy - I_zz) * q * r / I_xx + tau_phi / I_xx
    q_dot = (I_zz - I_xx) * p * r / I_yy + tau_theta / I_yy
    r_dot = (I_xx - I_yy) * p * q / I_zz + tau_psi / I_zz

    v_dot = np.array([[p_dot], [q_dot], [r_dot]])
    
    return v_dot

def get_n_dot_dot(X: np.ndarray, u: np.ndarray, I: np.ndarray, lift_constant: float, 
        motor_torque_constant: float, arm_length: float) -> np.ndarray:
# X = [[x], [y], [z], [x_dot], [y_dot], [z_dot], [phi], [theta], [psi], [phi_dot], [theta_dot], [psi_dot]]
# u = [[T], [t_phi], [t_theta], [t_psi]]
# I = [[I_xx, 0, 0], [0, I_yy, 0], [0, 0, I_zz]]
    phi = X[6, 0]
    theta = X[7, 0]
    psi = X[8, 0]
    phi_dot = X[9, 0]
    theta_dot = X[10, 0]
    psi_dot = X[11, 0]
    
    I_xx = I[0, 0]
    I_yy = I[1, 1]
    I_zz = I[2, 2]
    
    v = get_v(phi, theta, psi, phi_dot, theta_dot, psi_dot)

    v_dot = get_v_dot(phi, theta, psi, phi_dot, theta_dot, psi_dot, I, u, lift_constant, 
                      motor_torque_constant, arm_length)
    
    W_minus_one_n = get_W_minus_one_n(phi, theta, psi)

    W_minus_one_n_derivative = get_W_minus_one_n_derivative(phi, theta, psi, phi_dot, theta_dot, psi_dot)
    
    # equation (12)
    n_dot_dot = W_minus_one_n_derivative @ v + W_minus_one_n @ v_dot

    return n_dot_dot

def get_epsilon_dot_dot(X: np.ndarray, u: np.ndarray, m: float) -> np.ndarray:
# X = [[x], [y], [z], [x_dot], [y_dot], [z_dot], [phi], [theta], [psi], [phi_dot], [theta_dot], [psi_dot]]
# u = [[T], [t_phi], [t_theta], [t_psi]]
# m = mass of quad
    g = -9.81
    T = u[0 ,0]

    phi = X[6, 0]
    theta = X[7, 0]
    psi = X[8, 0]
    
    # equation (10)
    x_dot_dot = (T / m) * (np.cos(psi) * np.sin(theta) * np.cos(phi) + np.sin(psi) * np.sin(phi))
    y_dot_dot = (T / m) * (np.sin(psi) * np.sin(theta) * np.cos(phi) - np.cos(psi) * np.sin(phi))
    z_dot_dot = g + (T / m) * np.cos(theta) * np.cos(phi)
    
    epsilon_dot_dot = np.array([[x_dot_dot], [y_dot_dot], [z_dot_dot]])
    return epsilon_dot_dot

def quad_model(t: np.float64, X: np.ndarray, u: np.ndarray) -> np.ndarray: 
# reference https://sal.aalto.fi/publications/pdf-files/eluu11_public.pdf
# input state: 
# X = [[x], [y], [z], [x_dot], [y_dot], [z_dot], [phi], [theta], [psi], [phi_dot], [theta_dot], [psi_dot]]
# u = [[T], [tau_phi], [tau_theta], [tau_psi]]

#output state: 
#X_dot = [[x_dot], [y_dot], [z_dot], [x_dot_dot], [y_dot_dot], [z_dot_dot], [phi_dot], [theta_dot], [psi_dot], [phi_dot_dot], [theta_dot_dot], [psi_dot_dot]]
    assert(X.shape == (12,1))
    
    mass = 1
    I = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    lift_constant = 1
    motor_torque_constant = 1
    arm_length = 1
    
    x = X[0, 0]
    y = X[1, 0]
    z = X[2, 0]
    x_dot = X[3, 0]
    y_dot = X[4, 0]
    z_dot = X[5, 0]
    phi = X[6, 0]
    theta = X[7, 0]
    psi = X[8, 0]
    phi_dot = X[9, 0]
    theta_dot = X[10, 0]
    psi_dot = X[11, 0]
    
    epsilon_dot_dot = get_epsilon_dot_dot(X, u, mass)
    n_dot_dot = get_n_dot_dot(X, u, I, lift_constant, motor_torque_constant, arm_length)
    
    X_dot = np.zeros((12, 1))
    X_dot[0, 0] = x_dot
    X_dot[1, 0] = y_dot
    X_dot[2, 0] = z_dot
    X_dot[3, 0] = epsilon_dot_dot[0, 0]
    X_dot[4, 0] = epsilon_dot_dot[1, 0]
    X_dot[5, 0] = epsilon_dot_dot[2, 0]
    X_dot[6, 0] = phi_dot
    X_dot[7, 0] = theta_dot
    X_dot[8, 0] = psi_dot
    X_dot[9, 0] = n_dot_dot[0, 0]
    X_dot[10, 0] = n_dot_dot[1, 0]
    X_dot[11, 0] = n_dot_dot[2, 0]

    return X_dot 

time_start = 0
time_end = 10

# NOTE initial state vector needs to be shape (12,) i.e. have only 1 dimension. 
# But the value that gets passed into quad_model will have shape (12,1)
y0 = np.array([0] * 12)

# NOTE: "u" can be a 4 x N matrix where column i is the control input at time step i and each column is: [[T], [tau_phi], [tau_theta], [tau_psi]]
u = np.array([[0], [0], [0], [0]]) 

# tuple of arguments that get passed to quad_model. Add arguments as needed
args = (u, 1)

sol = solve_ivp(quad_model, [time_start, time_end], y0, vectorized = True, dense_output = True, args = args)

plot_results(sol)
