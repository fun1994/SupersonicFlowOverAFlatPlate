# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 11:41:45 2024

@author: HFKJ059
"""

import numpy as np
from matplotlib import pyplot as plt


def read_1d(path, filename):
    with open("./data/" + path + "/" + filename + ".txt", "r") as file:
        data = file.read()
    data = data.split()
    for i in range(len(data)):
        data[i] = float(data[i])
    data = np.array(data)
    return data

def read_3d(path, filename, len_x, len_y):
    data = []
    with open("./data/" + path + "/" + filename + ".txt", "r") as file:
        while True:
            line = file.readline()
            if not line:
                break
            data_temp = line.split()
            data.append([])
            for i in range(len_x):
                data[-1].append(data_temp[i*len_y:(i+1)*len_y])
    for i in range(len(data)):
        for j in range(len(data[i])):
            for k in range(len(data[i][j])):
                data[i][j][k] = float(data[i][j][k])
    data = np.array(data)
    return data

def read(path):
    x = read_1d(path, "x")
    y = read_1d(path, "y")
    t = read_1d(path, "time")
    rho = read_3d(path, "rho", len(x), len(y))
    u = read_3d(path, "u", len(x), len(y))
    v = read_3d(path, "v", len(x), len(y))
    T = read_3d(path, "T", len(x), len(y))
    p = read_3d(path, "p", len(x), len(y))
    Ma = read_3d(path, "Ma", len(x), len(y))
    return x, y, t, rho, u, v, T, p, Ma

def plot(X, Y, data, title):
    plt.contourf(X, Y, data.T, levels=1000, cmap="jet")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(title)
    plt.colorbar()
    plt.show()


x1, y1, t1, rho1, u1, v1, T1, p1, Ma1 = read("isothermal")
X1, Y1 = np.meshgrid(x1, y1)
plot(X1, Y1, rho1[-1, :, :], "rho under isothermal wall")
plot(X1, Y1, T1[-1, :, :], "T under isothermal wall")
plot(X1, Y1, p1[-1, :, :], "p under isothermal wall")
plot(X1, Y1, Ma1[-1, :, :], "Ma under isothermal wall")

x2, y2, t2, rho2, u2, v2, T2, p2, Ma2 = read("adiabatic")
X2, Y2 = np.meshgrid(x2, y2)
plot(X2, Y2, rho2[-1, :, :], "rho under adiabatic wall")
plot(X2, Y2, T2[-1, :, :], "T under adiabatic wall")
plot(X2, Y2, p2[-1, :, :], "p under adiabatic wall")
plot(X2, Y2, Ma2[-1, :, :], "Ma under adiabatic wall")
