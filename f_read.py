from scipy import *

def f_data(file,nx,ny):
    F0 = zeros((nx,ny))
    F1 = zeros((nx,ny))
    F2 = zeros((nx,ny))
    W = zeros((nx,ny))
    Z = zeros((nx,ny))
    i=0;j=0
    for line in open(file):
        if len(line)>1:
            line = line.split()
            F1[i,j] = line[2]
            F2[i,j] = line[3]
            F0[i,j] = line[4]
            Z[i,j] = line[5]
            W[i,j] = line[6]
            i+=1
        else:
            j+=1
            i=0
    return F1, F2, F0, Z, W

def f_fixed_r_data(file,ny):
    F0 = zeros(ny)
    F1 = zeros(ny)
    F2 = zeros(ny)
    W = zeros(ny)
    Z = zeros(ny)
    i=0
    for line in open(file):
        if len(line)>1:
            line = line.split()
            F1[i] = line[1]
            F2[i] = line[2]
            F0[i] = line[3]
            Z[i] = line[4]
            W[i] = line[5]
            i+=1
    return F1, F2, F0, Z, W

def ff_data(file,nx,ny):
    file = open(file)
    Ttot = zeros((nx,ny))
    rho = zeros((nx,ny))
    T34 = zeros((nx,ny))
    for j in range(ny):
        for i in range(nx):
            line = file.readline().split()
            Ttot[i,j] = line[0]
            rho[i,j] = line[1]
            T34[i,j] = line[2]
    return Ttot, rho, T34
