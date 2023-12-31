#https://github.com/yotamashkenazy/try.git
def optimized_derivative(t, states_n, planets, r_max, m_p, m_b=1.98e30):
    G = 6.6743e-11  # gravitational constant [m^3 kg^-1 s^-1]
    c = 299792458   # speed of light [m/s]

    x_ab = -states_n[:3]
    v_a = states_n[3:]
    
    r_ab = (x_ab@x_ab)**0.5
    
    
    dis = planets + x_ab
    dis_mag = np.linalg.norm(dis, axis=1)
    cond = dis_mag < r_max
    a_pl = np.sum(G * m_p * dis[cond] * (dis_mag[cond] ** -3)[:, np.newaxis], axis=0)
    
    F_g = x_ab*G * m_b / r_ab**3
    a_pn = F_g *(v_a@v_a - 4 * G * m_b / r_ab) - 4 * v_a * v_a@F_g
    return np.append(v_a,F_g + a_pn/c**2 + a_pl)
# Define initial conditions
# Format: [x, y, z, vx, vy, vz]
initial_conditions = np.array([1.0e9, 0.0, 0.0, 0.0, 500.0, 0.0])  # Example values in meters and meters per second

# Define planetary data with more planets
# Example: Earth and additional planets
planets = np.array([
    [0.0, 0.0, 0.0],        # Sun's position (centered at the origin)
    [1.5e11, 0.0, 0.0],     # Earth's position (approximately 1 AU from the Sun)
    [2.5e11, 0.0, 0.0],     # Planet 2's position
    [3.5e11, 0.0, 0.0],     # Planet 3's position
    # Add more planets as needed
])

# Other parameters
r_max = 1.0e11  # Maximum distance for calculations (arbitrary)

# Masses of the bodies
m_p = 5.972e24  # Earth's mass in kg
m_b = 1.989e30  # Sun's mass in kg

for _ in range(100000):
    optimized_derivative(0, initial_conditions, planets, r_max, m_p, m_b)
