import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import OdeSolution




def plot_results(sol: OdeSolution) -> None:
    # TODO plot results from sol
    return

def get_W_minus_one_n_derivative(phi: float, theta: float, psi: float) -> np.ndarray:
    # equation (12)
    zero_zero = 0
    zero_one = phi_dot * np.cos(phi) * np.tan(theta) + theta_dot * np.sin(phi) / (np.cos(theta)**2)
    zero_two = -phi_dot * np.sin(phi) * np.cos(theta) + 
    one_zero = 0
    one_one = 
    one_two = 
    two_zero = 0
    two_one = 
    two_two = 
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
    v = W_n @ (np.array([phi_dot, theta_dot, psi_dot]).reshape((3, 1)))
    return v

def get_W_minus_one_n(phi:float, theta:float, psi:float) -> np.ndarray:
    # equation (4)
    W_minus_one_n = np.array([[1, np.sin(phi) * np.tan(theta), np.cos(phi) * np.tan(theta)],
                              [0, np.cos(phi), -np.sin(phi)],
                              [0, np.sin(phi) / np.cos(theta), np.cos(phi) / np.cos(theta)]])
    return W_minus_one_n

def get_u_from_omega(omega: np.ndarray, lift_constant: float, motor_torque_constant: float, 
                     arm_length: float) -> np.ndarray:
    # NOTE: counter-clockwise torque is positive, clockwise torque is negative
    # equations (7) and (8) # 
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

def get_omega_tau(u: np.ndarray, lift_constant: float, motor_torque_constant: float, arm_length: float) -> float:
    T = u[0, 0]
    tau_phi = u[1, 0]
    tau_theta = u[2, 0]
    tau_psi = u[3, 0]
    
    omega = get_omega_from_u(u)
    return np.sum(omega)

def get_v_dot(phi: float, theta: float, psi:: float, phi_dot: float, theta_dot: float, 
        psi_dot: float, I: np.ndarray, u: np.ndarray, lift_constant: float) -> np.ndarray:
    # equation (11)
    I_xx = I[0, 0]
    I_yy = I[1, 0]
    I_zz = I[2, 0]

    tau_phi = u[1]
    tau_theta = u[2]
    tau_psi = u[2]
    
    v = get_v(phi, theta, psi, phi_dot, theta_dot, psi_dot)
    p = v[0, 0]
    q = v[1, 0]
    r = v[2, 0]
    
    omega_tau = get_omega_tau(u)

    p_dot = (I_yy - I_zz) * q * r / I_xx - I_xx * (q / I_xx) * omega_tau + tau_phi / I_xx
    q_dot = (I_zz - I_xx) * p * r / I_yy - I_yy * (-p / I_yy) * omega_tau + tau_theta / I_yy
    r_dot = (I_xx - I_yy) * p * q / I_zz - I_zz * 0 * omega_tau + tau_psi / I_zz

    v_dot = np.array([[p_dot], [q_dot], [r_dot]])
    return v_dot

def get_n_dot_dot(X: np.ndarray, u: np.ndarray, I: np.ndarray) -> np.ndarray:
# X = [x, y, z, x_dot, y_dot, z_dot, phi, theta, psi, phi_dot, theta_dot, psi_dot]   
# u = [T, t_phi, t_theta, t_psi]
# I = [[I_xx, 0, 0], [0 I_yy 0], [0 0 I_zz]]
    phi = X[6]
    theta = X[7]
    psi = X[8]
    phi_dot = X[9]
    theta_dot = X[10]
    psi_dot = X[11]
    
    I_xx = I[0, 0]
    I_yy = I[1, 1]
    I_zz = I[2, 2]
    
    v = get_v(phi, theta, psi, phi_dot, theta_dot, psi_dot)

    v_dot = get_v_dot(phi, theta, psi, phi_dot, theta_dot, psi_dot, I)
    
    W_minus_one_n = get_W_minus_one_n(phi, theta, psi)

    W_minus_one_n_derivative = get_W_minus_one_n_derivative(phi, theta, psi)
    
    n_dot_dot = W_minus_one_n_derivative @ v + W_minus_one_n @ v_dot

    return n_dot_dot

def get_epsilon_dot_dot(X: np.ndarray, u: np.ndarray, m: float) -> np.ndarray:
# X = [x, y, z, x_dot, y_dot, z_dot, phi, theta, psi, phi_dot, theta_dot, psi_dot]
# u = [T, t_phi, t_theta, t_psi]
# m = mass of quad
    g = -9.81
    T = u[0 ,0]

    phi = X[6]
    theta = X[7]
    psi = X[8]
    
    # equation (10)
    x_dot_dot = (T / m) * (np.cos(psi) * np.sin(theta) * np.cos(phi) + np.sin(psi) * np.sin(phi))
    y_dot_dot = (T / m) * (np.sin(psi) * np.sin(theta) * np.cos(phi) - np.cos(psi) * np.sin(phi))
    z_dot_dot = g + (T / m) * np.cos(theta) * np.cos(phi)
    
    epsilon_dot_dot = np.array([x_dot_dot, y_dot_dot, z_dot_dot])
    return epsilon_dot_dot

def quad_model(t: np.ndarray, X: np.ndarray): # TODO add control input "u" as parameter to quad_model
# reference https://sal.aalto.fi/publications/pdf-files/eluu11_public.pdf
# input state: 
# X = [x, y, z, x_dot, y_dot, z_dot, phi, theta, psi, phi_dot, theta_dot, psi_dot]

#output state: 
#X_dot = [x_dot, y_dot, z_dot, x_dot_dot, y_dot_dot, z_dot_dot, phi_dot, theta_dot, psi_dot, phi_dot_dot, theta_dot_dot, psi_dot_dot]
    m = 1
    I = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

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
    
    epsilon_dot_dot = get_epsilon_dot_dot(X, u, m)
    n_dot_dot = get_n_dot_dot(X, u, I)
    
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
