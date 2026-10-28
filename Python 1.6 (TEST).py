import numpy as np
import matplotlib.pyplot as plt


#Parameters
# Geometry & Printing
H = 0.0004               # Filament height and width (m)
L_total = 0.1            # Total length of the ellipse (m)
v_print = 0.03           # Printing speed (m/s) -> 30 mm/s

# Thermal properties (PLA)
rho = 1240.0             # Density (kg/m^3)
cp = 1800.0              # Heat capacity (J/kg*K)
k = 0.13                 # Thermal conductivity (W/m*K)
alpha = k/(rho*cp)   # Thermal diffusivity (m^2/s)

# Temperatures (°C)
T_nozzle = 210.0
T_plate = 60.0
T_env = 20.0

# Heat transfer coefficients
h_plate = 150.0          # To the heated base plate (W/m^2*K)
h_air = 15.0             # To ambient air across 3 sides (W/m^2*K)
c1 = h_plate/(rho*cp*H)
c2 = 3*h_air/(rho*cp*H)


#Time
N = 100                  # Number of elements dividing the circle
ds = L_total/N        # Length of 1 element (m)


dt = ds/v_print        # Time to print 1 element (s)

t_print = N*dt         # Total printing time (s)
t_extra = 5.0            # Extra cooling time after printing (s)
steps_extra = int(t_extra / dt)

#Starting conditions and stopping conditions
n = 0           #Starting element
t = 0.0          #Starting time (s)
Biggest_difference = 5 #(°C)
T_list = [T_env]*(N+2)  # Starting temperatures list (°C)
T_dif = max(T_list) - min(T_list) #Starting difference (°C)

def update_Tempratures():
    for i in range(n):
        T_list[i] = T_list[i] + (c1*(T_plate - T_list[i]) + c2*(T_env - T_list[i]))*dt


while (n < N and T_dif > Biggest_difference):
    T_list_old = T_list.copy()
    if (n < N):
                n += 1
            else 
                n = n
    t += dt
    T_dif = max(T_list) - min(T_list)
    if (n < N):
         T_list_old[n+1] = T_nozzle
    else  
        T_list_old[n+1] = T_list_old[0]   

    for i in range(n):
            if (n < N):
                T_list_old[i+1] = T_nozzle
            else  
                T_list_old[i+1] = T_list_old[0]  
            T_list[i] = T_list_old[i] + (alpha*((T_list_old[i-1] - T_list_old[i] + T_list_old[i+1])/(ds**2)) + c1*(T_plate - T_list_old[i]) + c2*(T_env - T_list_old[i]))*dt
