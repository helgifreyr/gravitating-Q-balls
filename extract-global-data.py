from scipy import *
from scipy.integrate import simps,trapz,quad,cumtrapz
from scipy.interpolate import interp1d
import scipy.odr.odrpack as odrpack
from pylab import *
from f_read import *

m,freq,alpha,c1,c2,c3 = genfromtxt('res.txt')

r = genfromtxt('gridx.dat')
theta = genfromtxt('gridy.dat')

nx = len(r)
ny = len(theta)

############################################################################
# calculate the mass and angular momentum from derivatives at infinity
############################################################################
infF1x, infF2x, infF0x, infZx, infWx = f_fixed_r_data('fx-inf.txt',ny)

infM = sum(infF0x)/ny
print 'Mass from infinity: %f' % infM

infF1xx, infF2xx, infF0xx, infZxx, infWxx = f_fixed_r_data('fxx-inf.txt',ny)

infJ = sum(infWxx)/ny
print 'Angular momentum from infinity: %f' % infJ

############################################################################
# calculate the mass from asymptotics
############################################################################
F1, F2, F0, Z, W = f_data('funct.dat',nx,ny)

minF0 = amin(F0)
F0H = F0[0,0]
F1H = F1[0,0]
maxW = amax(W)
maxZ = amax(Z)

def f(B, x):
    return B[0]/x + B[1]/x**2

def w(B, x):
    return B[0]/x**2 + B[1]/x**3

Fmodel = odrpack.Model(f)
Wmodel = odrpack.Model(w)
sx=std(r[1:])
F0fits = zeros((ny,2))
F1fits = zeros((ny,2))
F2fits = zeros((ny,2))
Wfits = zeros((ny,2))

def fit(x,y,model,beta0=[1.0,1.0]):
    sx = std(x)
    sy = std(y)
    data = odrpack.RealData(x, y, sx=sx, sy=sy)
    odr = odrpack.ODR(data, model, beta0=beta0)
    output = odr.run()
    return output.beta

for i in range(ny):
    if i!=0:
        F0fits[i,:] = fit(r[1:],F0[1:,i],Fmodel,beta0=F0fits[i-1,:])
    else:
        F0fits[i,:] = fit(r[1:],F0[1:,i],Fmodel,beta0=[-infM,infJ])
    F1fits[i,:] = fit(r[1:],F1[1:,i],Fmodel)
    F2fits[i,:] = fit(r[1:],F2[1:,i],Fmodel)
    Wfits[i,:] = fit(r[1:],W[1:,i],Wmodel)

# plot(r[1:],F0[1:,-1],'k')
# plot(r[1:],f(F0fits[-1],r[1:]),'b')
# xlim(0,120)
# ylim(min(F0[1:,-1]),max(F0[1:,-1]))
# show()

# print F0fits[:,0]

asyM = -sum(F0fits,axis=0)[0]/ny
print 'Mass from asymptotics: %f' % asyM

############################################################################
# calculate the mass the stress-energy tensor: Smarr relations
############################################################################

Ttot, rho, T34 = ff_data('T44.dat',nx,ny)

Mint_r_theta = zeros((nx,ny))
Jint_r_theta = zeros((nx,ny))
Mint_theta = zeros(ny)
Mint_theta2 = zeros(ny)
Jint_theta = zeros(ny)
for j in range(1,ny-1):
    for i in range(1,nx-1):
        Mint_r_theta[i,j] = exp(F0[i,j] + 2*F1[i,j] + F2[i,j]) * r[i]**2 * Ttot[i,j]
        Jint_r_theta[i,j] = exp(F0[i,j] + 2*F1[i,j] + F2[i,j]) * r[i]**2 * T34[i,j]
    Mint_theta[j] = trapz(Mint_r_theta[:,j],r)
    # Mint_theta2[j] = quad(interp1d(r[1:-2],Mint_r_theta[1:-2,j]),min(r[1:-2]),max(r[1:-2]))[0]
    Jint_theta[j] = trapz(Jint_r_theta[:,j],r)

Mint_theta[0] = Mint_theta[1]
# Mint_theta2[0] = Mint_theta2[1]
Jint_theta[0] = Jint_theta[1]
Mint_theta[-1] = Mint_theta[-2]
# Mint_theta2[-1] = Mint_theta2[-2]
Jint_theta[-1] = Jint_theta[-2]
Mint = trapz(sin(theta)*Mint_theta,theta)
# Mint2 = quad(interp1d(theta,Mint_theta),min(theta),max(theta))[0]
Jint = trapz(sin(theta)*Jint_theta,theta)
print 'Mass from stress-energy tensor: %f' % Mint
# print 'Mass from stress-energy tensor: %f' % Mint2
print 'Angular momentum from stress-energy tensor: %f' % Jint


Zm, rm = genfromtxt('sup.dat')
output = [alpha, freq, c1, c2, c3, infM, infJ, Mint, Jint, minF0, F0H, F1H, maxW,
        Zm, rm]
tmp = open('tmp-py.txt','w')
for par in output:
    tmp.write(str(par)+' ')
