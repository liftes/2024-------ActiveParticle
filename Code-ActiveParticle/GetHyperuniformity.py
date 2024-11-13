import numpy as np
from numba import jit, njit, prange
from tqdm import trange

import norfetools as nt
from norfetools import plt

def generate_particles(Lx, Ly, r_c, a_plane_2, fraction, layout_type=1):
    effective_Ly = Ly * fraction
    a_plane_3 = r_c
    num_particles = 0

    if layout_type == 1:
        dx = a_plane_2
        dy = a_plane_2
        rows = int(np.floor(effective_Ly / dy))
        nx = int(np.floor((Lx / 2 - r_c) / dx) * 2 + 1)
        num_particles = nx * rows

        x_positions = np.random.uniform(r_c, Lx - r_c, num_particles)
        y_positions = np.random.uniform(0, effective_Ly, num_particles)
        particles = np.column_stack((x_positions, y_positions))

        # # 向左平移Lx/2
        # particles[:, 0] -= Lx / 2
        
        # # 对于小于0的横坐标加上Lx
        # particles[:, 0] = np.where(particles[:, 0] < 0, particles[:, 0] + Lx, particles[:, 0])

        return particles
    else:
        raise ValueError("Unsupported layout type. Only layout_type=1 (Square grid) is supported.")
@njit
def wca_potential(r, epsilon=1.0, sigma=1.0):
    if r < sigma:
        return 0.5 * (sigma - r)**2
    else:
        return 0.0

@njit
def compute_distance_matrix(particles, Lx, Ly):
    num_particles = len(particles)
    dist_matrix = np.zeros((num_particles, num_particles))
    for i in prange(num_particles):
        for j in prange(i+1, num_particles):
            dx = particles[i, 0] - particles[j, 0]
            dy = particles[i, 1] - particles[j, 1]
            # Apply periodic boundary conditions
            dx -= Lx * np.round(dx / Lx)
            dy -= Ly * np.round(dy / Ly)
            dist_matrix[i, j] = dist_matrix[j, i] = np.sqrt(dx**2 + dy**2)
    return dist_matrix

@njit
def total_energy_from_dist_matrix(dist_matrix, epsilon=1.0, sigma=1.0):
    energy = 0.0
    num_particles = len(dist_matrix)
    for i in prange(num_particles):
        for j in prange(i+1, num_particles):
            energy += wca_potential(dist_matrix[i, j], epsilon, sigma)
    return energy

@njit
def update_distance_matrix(dist_matrix, particles, particle_idx, new_position, Lx, Ly):
    num_particles = len(particles)
    for j in prange(num_particles):
        if j != particle_idx:
            dx = new_position[0] - particles[j, 0]
            dy = new_position[1] - particles[j, 1]
            # Apply periodic boundary conditions
            dx -= Lx * np.round(dx / Lx)
            dy -= Ly * np.round(dy / Ly)
            dist_matrix[particle_idx, j] = dist_matrix[j, particle_idx] = np.sqrt(dx**2 + dy**2)
    return dist_matrix

def monte_carlo_simulation(particles, Lx, Ly, r0, num_steps=100, epsilon=1.0, sigma=1.0, delta=0.1):
    sigma = r0
    num_particles = len(particles)
    dist_matrix = compute_distance_matrix(particles, Lx, Ly)
    current_energy = total_energy_from_dist_matrix(dist_matrix, epsilon, sigma)
    energyList = []
    # for _ in trange(num_steps, desc="MC Simulation", ncols=80):
    for _ in range(num_steps):
        for _ in range(particles.shape[0]):
            # Randomly select a particle
            particle_idx = np.random.randint(num_particles)
            old_position = particles[particle_idx].copy()
            
            # Propose a new position with a small displacement
            displacement = np.random.uniform(-delta, delta, 2)
            new_position = old_position + displacement
            
            # Apply periodic boundary conditions
            new_position[0] = new_position[0] % Lx
            new_position[1] = new_position[1] % Ly
            
            # Temporarily update distance matrix for the proposed move
            new_dist_matrix = dist_matrix.copy()
            new_dist_matrix = update_distance_matrix(new_dist_matrix, particles, particle_idx, new_position, Lx, Ly)
            
            # Calculate the energy change
            new_energy = total_energy_from_dist_matrix(new_dist_matrix, epsilon, sigma)
            delta_energy = new_energy - current_energy
            
            # Apply the Metropolis criterion
            if delta_energy < 0 or np.random.rand() < np.exp(-delta_energy):
                particles[particle_idx] = new_position
                dist_matrix = new_dist_matrix
                current_energy = new_energy
    
        print(f"{current_energy:.4f}", end=", ")
        energyList.append(current_energy)

    return particles, energyList

def save_particles_to_txt(particles, filename, Lx):
    # 向左平移Lx/2
    particles[:, 0] -= Lx / 2
    
    # 对于小于0的横坐标加上Lx
    particles[:, 0] = np.where(particles[:, 0] < 0, particles[:, 0] + Lx, particles[:, 0])

    with open(filename, 'w') as file:
        for particle in particles:
            line = f"{particle[0]} {particle[1]}\n"
            file.write(line)

# 示例调用
Lx = 50.0
Ly = 200.0
r_c = 0
a_plane_2 = 3.0
fraction = 0.8

particles = generate_particles(Lx, Ly, r_c, a_plane_2, fraction)

for r0 in np.arange(2.0, 4.2, 0.2):
    particles, energyList = monte_carlo_simulation(particles, Lx, Ly, r0, num_steps = 30)

    save_particles_to_txt(particles, f"Hyperuniformity-{r0:.2f}.txt", Lx)

    plt.figure(figsize=(3, 5))
    plt.scatter(particles[:, 0], particles[:, 1], s=1)
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    plt.xlim(0, Lx)
    plt.ylim(0, Ly)
    plt.gca().set_aspect('equal', adjustable='box')  # 确保比例符合真实比例
    nt.SaveFig(1, f"Loc-{r0:.2f}.png", "Figure/Init-Hyperuniformity/")

    plt.plot(energyList)
    nt.SaveFig(1, f"Energy-{r0:.2f}.png", "Figure/Init-Hyperuniformity/")
