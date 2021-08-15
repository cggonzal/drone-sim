from sympy import symbols, init_printing, Matrix, sin, cos, tan

init_printing()

# see here for model: https://www.cggonzalez.com/blog/model.html

x_dot, y_dot, z_dot = symbols("x_dot y_dot z_dot")
phi, theta, psi = symbols("phi theta psi")
phi_dot, theta_dot, psi_dot = symbols("phi_dot theta_dot psi_dot")
p, q, r = symbols("p q r")
tau_phi, tau_theta, tau_psi, f_t = symbols("tau_phi tau_theta tau_psi f_t")
m, g = symbols("m g")
I_x, I_y, I_z = symbols("I_x I_y I_z")
nu, nu_dot = symbols("nu nu_dot")

x_dot_dot = f_t * (cos(psi)*sin(theta)*cos(phi) + sin(psi)*sin(phi)) / m
y_dot_dot = f_t * (sin(psi)*sin(theta)*cos(phi) - cos(psi)*sin(phi)) / m
z_dot_dot = (f_t * cos(theta)*cos(phi) / m) - g

nu = Matrix([[p], [q], [r]])

nu_dot = Matrix([[(I_y - I_z) * q * r / I_x + tau_phi / I_x],
                 [(I_z - I_x) * p * r / I_y + tau_theta / I_y],
                 [(I_x - I_y) * p * q / I_z + tau_psi / I_z]])

nu_transformation_zero_one = phi_dot*cos(phi)*tan(theta) + theta_dot * sin(phi) / cos(theta)**2
nu_transformation_zero_two = -1 * phi_dot * sin(phi) * cos(theta) + theta_dot * cos(phi) / cos(theta)**2
nu_transformation_one_one = -1 * phi_dot * sin(phi)
nu_transformation_one_two = -1 * phi_dot * cos(phi)
nu_transformation_two_one = phi_dot * cos(phi) / cos(theta) + phi_dot * sin(phi) * tan(theta) / cos(theta)
nu_transformation_two_two = -1 * phi_dot * sin(phi) / cos(theta) + theta_dot * cos(phi) * tan(theta) / cos(theta)
nu_transformation = Matrix([[0, nu_transformation_zero_one, nu_transformation_zero_two],
                            [0, nu_transformation_one_one, nu_transformation_one_two],
                            [0, nu_transformation_two_one, nu_transformation_two_two]])

body_to_world_transformation = Matrix([[1, sin(phi) * tan(theta), cos(phi) * tan(theta)],
                                       [0, cos(phi), -1*sin(phi)],
                                       [0, sin(phi) / cos(theta), cos(phi) / cos(theta)]])

attitude = nu_transformation * nu + body_to_world_transformation * nu_dot
phi_dot_dot, theta_dot_dot, psi_dot_dot = attitude[0, 0], attitude[1, 0], attitude[2, 0]

# final result
X_dot = Matrix([[x_dot],
                [y_dot],
                [z_dot],
                [phi_dot],
                [theta_dot],
                [psi_dot],
                [x_dot_dot],
                [y_dot_dot],
                [z_dot_dot],
                [phi_dot_dot],
                [theta_dot_dot],
                [psi_dot_dot]])

print(X_dot[9, 0])
