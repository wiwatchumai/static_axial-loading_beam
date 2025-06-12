import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

t = sp.symbols('t')
x = sp.symbols('x')
u = sp.Function('u')
E, rho, A, l = sp.symbols('E rho A l')
PDE = sp.Eq(rho*u(x, t).diff(t, 2) - E*A*u(x, t).diff(x, 2), 0)

E = float(input("Enter the Young's modulus (E [Pascarl]): "))
rho = float(input("Enter the density (rho [kg/m^3]): "))
A = float(input("Enter the cross-sectional area (A [m^2]): "))
l = float(input("Enter the length of the beam (l [m]): "))

k = E * A / l

BC1 = sp.Eq(u(0, t), 0)
BC2 = sp.Eq(u(l, t).diff(x, 1), 0)
BC3 = sp.Eq(u(x, 0).diff(t, 1), 0)

i = int(input("Select the number of nodes : "))

forces = []
for n in range(1, i+1):
    F_n = float(input(f"Input force magnitude at node {n}: "))
    forces.append(F_n)

for n in range(1, i+1):
    u_n = sp.symbols(f'u{n}')
    u_np1 = sp.symbols(f'u{n+1}')
    F_n = k * (u_n - u_np1)
    print(f"F({n}) =")
    sp.pprint(F_n)

K = np.zeros((i, i))
for n in range(i - 1):
    K[n, n]     += k
    K[n, n+1]   -= k
    K[n+1, n]   -= k
    K[n+1, n+1] += k

F = np.array(forces)

K[0, :] = 0
K[0, 0] = 1
F[0] = 0

u = np.linalg.solve(K, F)

print("Nodal displacements (u_n):")
for n in range(i):
    print(f"u{n} = {u[n]} m.")

plt.figure()
plt.plot(range(i), u, marker='o')
plt.title("Nodal Displacement vs Node Number")
plt.xlabel("Node Number")
plt.ylabel("Displacement (m)")
plt.grid(True)
plt.show()
