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


def gamma(i):
    if i == 1:
        return np.kron(pauli(3), pauli(1))
    elif i == 2:
        return np.kron(pauli(3), pauli(2))
    elif i == 3:
        return np.kron(pauli(3), pauli(3))
    elif i == 4:
        return np.kron(pauli(1), np.eye(2))
    elif i == 5:
        return np.kron(pauli(2), np.eye(2))
    else:
        raise ValueError("Invalid input!")


def get_Hamiltonian(k1, k2, k3, k4, k5, u, b):
    H = np.sin(k1) * gamma(1) + np.sin(k2) * gamma(2) \
        + np.sin(k3) * gamma(3) + np.sin(k4) * gamma(4) \
        + (4 + u - np.cos(k1) - np.cos(k2) - np.cos(k3) - np.cos(k4) - np.cos(k5)) * gamma(5) \
        + (1j*b)/2 * (gamma(3) @ gamma(4) - gamma(4) @ gamma(3))
    return H


M = 3
u = 0.0
b = 0.5

varphi_2 = 0
varphi_3 = 0
varphi_4 = 0
varphi_5 = (1/2) * np.pi

k1_list = np.linspace(-np.pi * 0.05, np.pi * 0.05, num=1001, endpoint=True)
E = np.zeros((4, len(k1_list)))

for k1_index in range(len(k1_list)):
    k1 = k1_list[k1_index]
    k2 = M * k1 + varphi_2
    k3 = M**2 * k1 + varphi_3
    k4 = M**3 * k1 + varphi_4
    k5 = M**4 * k1 + varphi_5

    H = get_Hamiltonian(k1, k2, k3, k4, k5, u, b)
    E[:, k1_index], _ = la.eigh(H)

plt.figure(figsize=(6, 8))
plt.plot(k1_list/np.pi, E[0, :], color='black', linewidth=2)
plt.plot(k1_list/np.pi, E[1, :], color='black', linewidth=2)
plt.plot(k1_list/np.pi, E[2, :], color='black', linewidth=2)
plt.plot(k1_list/np.pi, E[3, :], color='black', linewidth=2)

plt.gca().set_xticklabels([])
plt.gca().set_yticklabels([])
for spine in plt.gca().spines.values():
    spine.set_linewidth(2)

plt.xlim(-0.05, 0.05)
plt.ylim(-6, 6)
plt.xticks((-0.05, 0, 0.05))
plt.yticks((-6, 0, 6))
plt.show()
