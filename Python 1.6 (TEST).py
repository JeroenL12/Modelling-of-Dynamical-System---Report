import numpy as np
import matplotlib.pyplot as plt


#Parameters
# Geometry & Printing
H = 0.0004               # Filament height and width (m)
L_total = 0.5            # Total length of the ellipse (m)
v_print = 0.03           # Printing speed (m/s)

# Thermal properties (PLA)
rho = 1240.0             # Density (kg/m^3)
cp = 1300.0              # Heat capacity (J/kg*K)
k = 0.13                 # Thermal conductivity (W/m*K)
alpha = k/(rho*cp)   # Thermal diffusivity (m^2/s)

# Temperatures (°C)
T_nozzle = 210.0
T_plate = 60.0
T_env = 20.0

# Heat transfer coefficients
h_plate = 150.0          # To the heated base plate (W/m^2*K)
h_air = 15.0             # To ambient air across 3 sides (W/m^2*K)
c1 = h_plate/(rho*cp*H)  # Coefficient 1
c2 = 3*h_air/(rho*cp*H)  # Coefficient 2


#Time
N = 1000                 # Number of elements dividing the circle
ds = L_total/N           # Length of 1 element (m)
dt = ds/v_print          # Time to print 1 element (s)           

#Starting conditions and stopping conditions
n = 0                    #Starting element
t = 0.0                  #Starting time (s)
T_dif_allowed = 0.5      #(°C)
T_list = [0]*(N)         # Starting temperatures list (°C)
T_dif = 0                #Starting difference (°C)
a = 0                    #state that print is at last element (Print finished)

#Function for updating tempratures along the whole printed length
def update_Tempratures():
    for i in range(limit):
        
        if (i == 0 and n < N and n != 1):
            T_left = T_list_old[i]
            T_right = T_list_old[i+1]

        elif (i == 0 and n == N):
            T_left = T_list_old[n-1]
            T_right = T_list_old[i+1]

        elif((i == n-1 and n < N) or n==1):
            T_list[i] = T_nozzle
            continue

        elif(i == n-1 and n == N):
            if a == 1:
                T_list[i] = T_nozzle
                continue
            else:
                T_right = T_list_old[0]
                T_left = T_list_old[i-1]

        else:
            T_left = T_list_old[i-1]
            T_right = T_list_old[i+1]
        
        T_list[i] = T_list_old[i] + (alpha*((T_left - 2*T_list_old[i] + T_right)/(ds**2)) + c1*(T_plate - T_list_old[i]) + c2*(T_env - T_list_old[i]))*dt

#Loop for finishing the whole print and checking if the difference is below a certain threshold
while (n < N or T_dif > T_dif_allowed):
    T_list_old = T_list.copy()
    
    if (n < N):
        n += 1
        limit = n
        if n == N:
            a=1
            
    else:
        a = 0
        limit = N
    
    t += dt
    if a == 1:
        print_time = t

    update_Tempratures()
    T_dif = max(T_list) - min(T_list)

    #If t becomes to big stop simulation
    if t > 10**4:
        print("Error: Maximum number of seconds exceeded 10000 (2.8 hours). Stopping simulation.")
        break

print("Total time: ", round(t, 4), "seconds")
print("Total printing time:", round(print_time, 4), "seconds")
print("Total cooling time: ", round(t-print_time, 4), "seconds")
