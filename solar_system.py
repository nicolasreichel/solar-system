import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import itertools as iter

# Constants
# G = 6.67408 * 10**(-11) # Gravitational constant in SI units
G = 39.478 # Gravitational constant in AU^3/(M_s * Yr^2)

# Parameters
dt = 0.002
tmax = 500
# tmax = 200
# stepsperframe = 1
numframes = int(tmax/dt)
# numframes = 100000
inter = 0
planet_list = []

class Planet(object):
    global N
    """A planet or other celestial object."""
    def __init__(self, name, color, mass, distance, pos_0, v_0):
        super(Planet, self).__init__()
        self.name = name
        self.color = color
        self.mass = mass # Mass in Solar Masses
        self.au = distance # Distance from the sun in AU
        self.pos = [np.array(pos_0)] # Position in xy-plane [[x,y],...] [AU]
        # self.pos = pos_0 # Initial position vector [x_0,y_0] [AU]
        self.v_0 = np.array(v_0) # Initial velocity vector [v_x,v_y] [AU/Yr]
        self.vel = [np.array(v_0)] # Initial velocity vector [v_x,v_y] [AU/Yr]

# Calculate gravitational force on planet i from planet j
def F_g(planet_i,planet_j):
    m_j = planet_j.mass
    m_i = planet_i.mass
    pos_j = planet_j.pos[-1]
    pos_i = planet_i.pos[-1]
    r = (pos_j-pos_i)
    # DEBUG: distances
    # print("distance = ",r)
    # print("norm = ",np.linalg.norm(r))
    if np.linalg.norm(r) == 0: # Edge case
        print("ERROR COLLISION")
        return 0
    else:
        return (G*m_j*m_i/(np.linalg.norm(r)**3) * r)

# Generate force matrix
def gen_f_ij(planet_list):
    # Create a matrix f_ij where f[i,j] is the force on planet i from planet j.
    # Total force on planet i is sum(f[i,:])
    f_ij = np.zeros((len(planet_list),len(planet_list)),dtype=object)
    for i in range(0,len(planet_list)):
        for j in range(i+1,len(planet_list)):
            planet_i = planet_list[i]
            planet_j = planet_list[j]
            f_ij[i,j] = F_g(planet_i,planet_j)
            f_ij[j,i] = -1 * F_g(planet_i,planet_j)
    # DEBUG: look at forces
    # print(f_ij[1,:])
    return f_ij

# Perform one step of verlet integration on position vector x. Global error order O(t*e^at * dt^2)
def verlet(planet_list,dt):
    f_ij = gen_f_ij(planet_list)
    # DEBUG: force matrix
    # print(f_ij)
    for i,planet in enumerate(planet_list):
        # Acceleration on planet i equals sum of all forces acting on it, divided by its mass
        a_g = sum(f_ij[i,:])/planet.mass
        # DEBUG: gravity force
        # print("a_g = ",a_g)
        # print("Verlet, pos = ",planet.pos)
        if len(pos) == 1:
            # DEBUG:
            # print("verlet 1")
            x_old = planet.pos[0]
            x_new = x_old + v_0*dt + 0.5*a_g*dt # local error O(dt^3)
        else:
            # DEBUG:
            # print("verlet 2")
            x_old = planet.pos[-1]
            x_older = planet.pos[-2]
            # DEBUG: func thingy
            # print(pos)
            # print(a_g)
            x_new = 2*x_old - x_older + a_g*(dt**2) # local error order O(dt^4)
        # Approximate velocity with forward difference
        v_new = (x_new-x_old)/dt # local error O(dt^2)
        planet.vel.append(v_new)
        # Add new position to vector
        planet.pos.append(x_new)
    return

def leapfrog(planet_list,dt):
    for i,planet in enumerate(planet_list):
        x_old = planet.pos[-1]
        v_old = planet.vel[-1]
        if len(planet.pos) == 1:
            x_new = x_old + 0.5*v_old*dt
        else:
            x_new = x_old + v_old*dt
        planet.pos.append(x_new)
    # Calculate new velocities by acceleration at new point
    f_ij = gen_f_ij(planet_list)
    for i,planet in enumerate(planet_list):
        v_old = planet.vel[-1]
        a_g = sum(f_ij[i,:])/planet.mass
        v_new = v_old + a_g*dt
        planet.vel.append(v_new)

def animate(step):
    global planet_list,dt
    # Do one step of integration
    # verlet(planet_list,dt)
    leapfrog(planet_list,dt)
    # Retrieve position data
    planet_pos_list = []
    planet_color_list = []
    for planet in planet_list:
        # Retrieve latest position and add to vector
        # DEBUG: sun flying away
        # print("pos = ", planet.pos[-1])
        planet_pos_list.append(planet.pos[-1])
        # planet_color_list.append(planet.color)
    # Update scatter plot
    # print("positions: ",planet_pos_list)
    scat.set_offsets(planet_pos_list)
    time_text.set_text('{time:0.4f} Yr, dt={dt}'.format(time=step*dt,dt=dt))
    return scat,time_text

# CREATE PLANETS
# Name
Sun =          Planet("Sun",     '#ffcd00',                     1,     0, [0,    0], [0,0])
Mercury =      Planet("Mercury", '#dfddd7',    0.16601 * 10**(-6),  0.39, [0, 0.39], [10.0,0])
# Mercury =      Planet("Mercury", '#dfddd7',  0.16601 * 10**(-6),  0.39, [0, 0.39], [0,0])
Venus =        Planet("Venus", '#e5c568',       2.4478 * 10**(-6), 0.723, [0,0.723], [7.380,0])
Earth =        Planet("Earth", '#0671cb',       3.0035 * 10**(-6),     1, [0,    1], [6.28,0])
Mars =         Planet("Mars", '#b9480c',       0.32271 * 10**(-6), 1.524, [0,1.524], [5.080,0])
Jupiter =      Planet("Jupiter", '#d6ce84',     954.79 * 10**(-6), 5.203, [0,5.203], [2.76,0])
Saturn =       Planet("Saturn", '#f7eb78',      285.88 * 10**(-6), 9.539, [0,9.539], [2.05,0])
Uranus =       Planet("Uranus", '#b8e9f5',      43.662 * 10**(-6), 19.18, [0,19.18], [1.43,0])
Neptune =      Planet("Neptune", '#68b8f7',     51.514 * 10**(-6), 30.06, [0,30.06], [1.14,0])
Cheeky_Pluto = Planet("Pluto", '#fdd689',     0.007396 * 10**(-6), 39.53, [0,39.53], [0.991,0])
Asteroid =     Planet("4 Vesta", '#9f4400', 1,  None, [40,  40], [-40,-39.95]) # 4 Vesta
# Black_Hole =     Planet('#9f4400',    100,  None, [40,  40], [-30,-20])

# planet_list = [Sun,Mercury,Venus,Earth,Mars,Jupiter,Saturn,Uranus,Neptune,Cheeky_Pluto]
planet_list = [Sun,Mercury,Venus,Earth,Mars,Jupiter,Saturn,Uranus,Neptune,Cheeky_Pluto, Asteroid]
# planet_list = [Sun,Jupiter]

# plt.grid()
fig = plt.figure()
ax = plt.axes(xlim=(-50, 50), ylim=(-50, 50))
# ax = plt.axes(xlim=(-5, 5), ylim=(-5, 5))
ax.set_facecolor('k')
plt.gca().set_aspect('equal', adjustable='box')
plt.tight_layout()

# Retrieve initial data
planet_pos_list = [[],[]]
planet_color_list = []
for planet in planet_list:
    planet_pos = planet.pos[0]
    planet_pos_list[0].append(planet_pos[0])
    planet_pos_list[1].append(planet_pos[1])
    planet_color_list.append(planet.color)
scat = ax.scatter(planet_pos_list[0],planet_pos_list[1],c=planet_color_list,s=4)
# DEBUG: initial pos
# print(planet_pos_list)
time_text = ax.text(-45,45,'0 Yr, dt={dt}'.format(dt=dt),color='white')

anim = animation.FuncAnimation(fig, animate, frames=numframes, interval=inter, blit=True, repeat=False)
# anim.save("solar_system_animation_t10.gif", writer='imagemagick',fps=30)
plt.show()
