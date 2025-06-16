import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

# Variables definittion for a beam in two-dimensional space: x-t (displacement-time)
t=sp.symbols('t')
x=sp.symbols('x')
# Define the displacement function of a beam
u=sp.Function('u')
E, rho, A, l = sp.symbols('E rho A l')  # Young's modulus, density, cross-sectional area, length
# Define the wave equation for a beam
PDE = sp.Eq(rho*u(x, t).diff(t, 2) - E*A*u(x, t).diff(x, 2), 0)

E=float(input("Enter the Young's modulus (E [Pascarl]): "))
rho=float(input("Enter the density (rho [kg/m^3]): "))
A=float(input("Enter the cross-sectional area (A [m^2]): "))
l=float(input("Enter the length of the beam (l [m]): "))

k = E * A / l  # Structural rigidity

# Define the BCs
BC1 = sp.Eq(u(0, t), 0)  # Boundary condition at x=0
BC2 = sp.Eq(u(l, t).diff(x, 1), 0)  # Boundary condition at x=l
BC3 = sp.Eq(u(x, 0).diff(t, 1), 0)  # Boundary condition at t=0

# Input the number of nodes: i
i = int(input("Select the number of nodes : "))

forces = []
for n in range(i):
    if n == i - 1:
        F_n = float(input(f"Input force magnitude at node {n+1}: "))
    else:
        F_n = 0.0
    forces.append(F_n)

for n in range(1, i+1):
    u_n = sp.symbols(f'u{n}')
    u_np1 = sp.symbols(f'u{n+1}')
    F_n = k * (u_n - u_np1)
    print(f"F({n}) =")
    sp.pprint(F_n)

# Assemble global stiffness matrix
K = np.zeros((i, i))
for n in range(i - 1):
    K[n, n]     += k
    K[n, n+1]   -= k
    K[n+1, n]   -= k
    K[n+1, n+1] += k

# Convert forces list to numpy array
F = np.array(forces)

# Apply boundary condition: u0 = 0 (fixed at node 0)
K[0, :] = 0
K[0, 0] = 1
F[0] = 0

# Solve for displacements
u = np.linalg.solve(K, F)

print("Nodal displacements (u_n):")
for n in range(i):
    print(f"u{n} = {u[n]} m.")

# Plot displacement vs node number
plt.figure()
plt.plot(range(i), u, marker='o')
plt.title("Nodal Displacement vs Node Number")
plt.xlabel("Node Number")
plt.ylabel("Displacement (m)")
plt.grid(True)
plt.show()
