# -*- coding: utf-8 -*-
import numpy as np
from numpy import linalg as la
import matplotlib.pyplot as plt


def pauli(i):
    if i == 0:
        return np.array([[1, 0], [0, 1]])
    elif i == 1:
        return np.array([[0, 1], [1, 0]])
    elif i == 2:
        return np.array([[0, -1j], [1j, 0]])
    elif i == 3:
        return np.array([[1, 0], [0, -1]])
    else:
        raise ValueError("Invalid input!")


def get_Hamiltonian(k1, k2, k3, u):
    H = np.sin(k1) * pauli(1) + \
        np.sin(k2) * pauli(2) + \
        (2 + u - np.cos(k1) - np.cos(k2) - np.cos(k3)) * pauli(3)
    return H


M = 3
u = 0.0
varphi_2 = 0.0 * np.pi
varphi_3 = (1/2) * np.pi

k1_list = np.linspace(-np.pi, np.pi, num=801, endpoint=True)
E = np.zeros((2, len(k1_list)))

for k1_index in range(len(k1_list)):
    k1 = k1_list[k1_index]
    k2 = M * k1 + varphi_2
    k3 = M**2 * k1 + varphi_3
    H = get_Hamiltonian(k1, k2, k3, u)
    E[:, k1_index], _ = la.eigh(H)

plt.figure(figsize=(6, 8))
plt.plot(k1_list/np.pi, E[0, :], color='black', linewidth=2)
plt.plot(k1_list/np.pi, E[1, :], color='black', linewidth=2)

plt.gca().set_xticklabels([])
plt.gca().set_yticklabels([])
for spine in plt.gca().spines.values():
    spine.set_linewidth(2)

plt.xlim(-1, 1)
plt.ylim(-6, 6)
plt.xticks((-1, 0, 1))
plt.yticks((-6, 0, 6))
plt.show()
