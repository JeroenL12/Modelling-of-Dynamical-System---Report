import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#Parameters

# Motor parameters
r_m = 0.0064       # Motor pulley radius [m]
L_m = 0.003      # Motor inductance [H]
R_m = 1.8        # Motor resistance [Ohm]
J_m = 5.7e-6       # Motor rotor inertia [kg*m^2]
K_m = 0.22       # Motor constant [Nm/A]

#3d printer parameters
m = 0.6          # Tool head mass [kg]
J_tot = J_m + m*(r_m**2)/2  # Total effective inertia of the system
Max_x = 0.22       # Maximum x position [m]
Max_y = 0.22       # Maximum y position [m]


#Friction parameters
mu_r = 0.003     # Bearing rotational friction [Nm*s/rad]
mu_lx = 3      # Linear friction x-axis [N*s/m]
mu_ly = 3      # Linear friction y-axis [N*s/m]

#Time simulation parameters
total_time = 1 #Total simulation time [s]
step_size = 0.01 #Time step size [s]

#Reference trajectory parameters
omega = 2*np.pi # Frequency of the reference trajectory [rad/s]
a = (2/6)*Max_x  # Amplitude of the reference trajectory [m]
b = (1/6)*Max_y  # Amplitude of the reference trajectory [m]
A = (J_tot)/(K_m*r_m)
B = mu_r/J_tot + (mu_lx*(r_m**2))/(2*J_tot)
C = mu_r/J_tot + (mu_ly*(r_m**2))/(2*J_tot)

#Input Voltage Function

def input_voltages(t):
    x_dot_ref = -a * omega * np.sin(omega * t)
    x_dubbeldot_ref = -a * (omega**2) * np.cos(omega * t)
    x_triple_dot_ref = a * (omega**3) * np.sin(omega * t)

    i_mx_ref = A*(x_dubbeldot_ref + B*x_dot_ref)
    i_mx_dot_ref = A*(x_triple_dot_ref + B*x_dubbeldot_ref)

    y_dot_ref = b * omega * np.cos(omega * t)
    y_dubbeldot_ref = -b * (omega**2) * np.sin(omega * t)
    y_triple_dot_ref = -b * (omega**3) * np.cos(omega * t)

    i_my_ref = A*(y_dubbeldot_ref + C*y_dot_ref)
    i_my_dot_ref = A*(y_triple_dot_ref + C*y_dubbeldot_ref)
    

    u_mx = L_m*i_mx_dot_ref + R_m*i_mx_ref + (K_m/r_m)*x_dot_ref
    u_my = L_m*i_my_dot_ref + R_m*i_my_ref + (K_m/r_m)*y_dot_ref

    return u_mx, u_my


#State-Space Function

def state_space(t, z):
    """
    State vector x:
      x[0] = x 
      x[1] = x_dot    
      x[2] = y   
      x[3] = y_dot   
      x[4] = i_mx       
      x[5] = i_my       
    """
    x, x_dot, y, y_dot, i_mx, i_my = z
    u_mx, u_my = input_voltages(t)
    
    #Electrical Dynamics
    di_mx = (u_mx - R_m * i_mx - (K_m/r_m)*x_dot) / L_m
    di_my = (u_my - R_m * i_my - (K_m/r_m)*y_dot) / L_m
    
    #Mechanical Accelerations

    x_dubbeldot = ((K_m*r_m)/(J_tot))*(i_mx) - ((mu_r/J_tot) + (mu_lx*(r_m**2))/(2*J_tot))*(x_dot)
    y_dubbeldot = ((K_m*r_m)/(J_tot))*(i_my) - ((mu_r/J_tot) + (mu_ly*(r_m**2))/(2*J_tot))*(y_dot)
    
    # Return dx/dt array
    return [
        x_dot,       
        x_dubbeldot, 
        y_dot,       
        y_dubbeldot, 
        di_mx,       
        di_my       
    ]

#Solve ODE's
x0 = [Max_x/2 + a, 0.0, Max_y/2, 0.0, 0.0, 0.0]
t_span = (0.0, total_time)
t_eval = np.linspace(t_span[0], t_span[1], int(total_time/step_size))

sol = solve_ivp(
    fun=state_space,
    t_span=t_span,
    y0=x0,
    t_eval=t_eval,
    method = 'Radau'
)

#Plots
plt.figure(figsize=(12, 10))

'''
# Motor Currents
plt.subplot(2, 2, 1)
plt.plot(sol.t, sol.y[4], label=r'$i_{m,x}$ (x current)')
plt.plot(sol.t, sol.y[5], label=r'$i_{m,y}$ (y current)')
plt.ylabel('Current [A]')
plt.grid(True)
plt.legend()

#Tool Head Position (x, y)
plt.subplot(2, 2, 2)
plt.plot(sol.t, sol.y[0], label='Tool Head X')
plt.plot(sol.t, sol.y[2], label='Tool Head Y')
plt.xlabel('Time [s]')
plt.ylabel('Position [m]')
plt.grid(True)
plt.legend()
'''
#Tool position (x, y)
ax3 =plt.subplot(1, 1, 1)
plt.plot(sol.y[0], sol.y[2], label='Tool Head Path')
plt.xlim(0, Max_x)
plt.ylim(0, Max_y)
ax3.set_box_aspect(1)
plt.xlabel('Position X [m]')
plt.ylabel('Position Y [m]')
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()