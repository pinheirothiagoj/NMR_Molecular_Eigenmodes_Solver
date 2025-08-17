###############################################################################
# RICE UNIVERSITY
# Department of Chemical and Biomolecular Engineering
# Author: Thiago Jose Pinheiro dos Santos, Ph.D.
###############################################################################

# This code computes the radial component of the diffusion propagator for two
# magnetic dipoles in real units. This is achieved by numerically solving
# the corresponding Fokker-Planck equation for the problem. The solution is
# then used to compute the multiexponential time decays (τ's) and their
# amplitudes P(τ)'s, as well as the dipole-dipole NMR autocorrelation function
# and 1/T1(0) = 1/T2(0). This code reads input data from 'input.txt'. If the
# code is executed in potential_def = "read", it also reads 'rdf.txt'.
#
# This code (A) can be used by itself for a single case evaluation or in
# conjunction with codes (B) and (C) for evaluations over different values of
# ri/rf or for single evaluations with NMR frequency dispersion, respectively.
#
# When using these codes or any of their parts, please cite:
#
# [1] T. J. Pinheiro dos Santos, B. Orcan-Ekmecki, W. G. Chapman, P. M. Singer,
# and D. N. Asthagiri, “A comprehensive molecular eigenmodes model for dipole-
# dipole NMR relaxation in real fluids”.
#
# [2] T. J. Pinheiro dos Santos, B. Orcan-Ekmecki, W. G. Chapman, P. M. Singer,
# and D. N. Asthagiri, “Theory and modeling of molecular modes in the NMR
# relaxation of fluids,” J. Chem. Phys., vol. 064108, p. 160, 2024.


import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.signal import savgol_filter
import sys

##############################################################################
# Setting up variables
##############################################################################

# Radial discretization

radial_disc = "edgebiased"  # "uniform" or "edgebiased"

# Potential of interaction between dipoles

potential_def = "read"  # "read" or "analytical"

# Input of other constants

filename = "input.txt"

def read_input_file(filename):
    with open(filename, 'r') as file:
        return [float(value) if '.' in value or 'e' in value or 'E' in value else int(value)
                for line in file
                if '=' in line and (value := line.split('=')[1].strip())]

T,D0,ri,rf,dr,gI,gS,S,Km,dx,Kp,req = read_input_file(filename)

##############################################################################
# Input of other parameters/constants
##############################################################################

n    = 2                         # NMR second moment
pi   = 3.14159265
kB   = 1.380649E-23              # units = J/K
mu0  = 1.2566370621219E-6        # units = N/A²
hbar = 6.62607015E-34/(2*np.pi)  # units = J.s

alpha = (1/4) * (mu0/(4*np.pi))**2 * hbar**2 * gI**2 * gS**2 * S * (S+1)  # units = m⁶/s²

# Radial discretization

Nr = len(np.arange(ri,rf+dr,dr))

# Root-finding range discretization

x = np.arange(5.0E7,9.0E12+dx, dx)  # units = √m
Nx = len(x)

# Time discretization

dt = 0.01e-12               # units = s
t  = np.arange(0,5e-10,dt)  # units = s
Nt = len(t)

##############################################################################
# Defining vectors and matrices
##############################################################################

A = np.zeros([Nr,Nr],float)
B = np.zeros([Nr],float)
Rs = np.zeros([Nr],float)
Bsol = np.zeros([Nr],float)

U = np.zeros([Nr],float)
dU = np.zeros([Nr],float)
d2U = np.zeros([Nr],float)

D = np.zeros([Nr],float)
dD = np.zeros([Nr],float)
d2D = np.zeros([Nr],float)
mu = np.zeros([Nr],float)
dmu = np.zeros([Nr],float)

Rsi = np.zeros([Nr],float)
Bsoli = np.zeros([Nr],float)

hi = np.zeros([Nx],float)
hf = np.zeros([Nx],float)

zr = np.zeros([Nr],float)
zx = np.zeros([Nx],float)

root = []
tau  = []
idx  = []
Rs   = []

Gt = np.zeros([Nt],float)

original_stdout = sys.stdout
np.set_printoptions(threshold=sys.maxsize)

##############################################################################
# Defining important functions and quantities
##############################################################################

def generate_grid_with_dr_min(ri, rf, Nr, dr_min, dr_center):
    
    L = rf - ri
    s = np.linspace(0, 1, Nr)  # normalized coordinate

    # Calculate a to satisfy dr(0) = dr(1) = dr_min
    a = (dr_min - dr_center) / (0 - 0.5)**2
    dr_s = a * (s - 0.5)**2 + dr_center
    
    # Normalize dr_s so that sum(dr_s) = L
    dr_s *= L / np.sum(dr_s)
    
    # Now integrate dr_s to get r grid points
    r = ri + np.cumsum(dr_s)
    r = np.insert(r, 0, ri)
    
    return r

# Radial parameters

if radial_disc == "uniform":

    r = np.arange(ri,rf+dr,dr)

elif radial_disc == "edgebiased":

    dr_min = dr / 100 # lowering 2 orders of magnitude at the edges
    r = generate_grid_with_dr_min(ri, rf, Nr-1, dr_min, dr)

else:

    print("ERROR: Unknown radial discretization flag: radial_disc.")
    sys.exit(1)

# Potential of interaction, diffusivity, and mobility

if potential_def == "analytical":
    
    for i in range(0,Nr):
        
        U[i] = Kp*(r[i]-req)**2
        dU[i] = 2*Kp*(r[i]-req)
        d2U[i] = 2*Kp
            
        D[i] = D0 * np.exp(-U[i]/(kB*T))
        dD[i] = -D0/(kB*T) * np.exp(-U[i]/(kB*T)) * dU[i]
        d2D[i] = (-D0/(kB*T)) * ( np.exp(-U[i]/(kB*T)) * d2U[i] + (-dU[i]/(kB*T)) * np.exp(-U[i]/(kB*T)) * dU[i] )

        mu[i] = D[i]/(kB*T)
        dmu[i] = dD[i]/(kB*T)

elif potential_def == "read":

    # Read file
    data = np.loadtxt("rdf.txt")
    r_vals = data[:, 0]    # unit = m
    rdf_vals = data[:, 1]  # dimensionless

    # Avoid zero or negative RDF values
    rdf_safe = np.clip(rdf_vals, 1e-12, None)

    # Compute potential of mean force (in J)
    U_vals = -np.log(rdf_safe) * kB * T

    # Create spline
    U_spline = CubicSpline(r_vals, U_vals, bc_type='natural', extrapolate=True)

    # Evaluate (vectorized)
    U = U_spline(r)

    dU = np.zeros_like(U)
    d2U = np.zeros_like(U)
    
    # Left boundary (forward difference)
    
    h = r[1] - r[0]
    
    dU[0] = (U[1] - U[0]) / h
    d2U[0] = (U[2] - 2*U[1] + U[0]) / (h**2)

    # Local spacing approximation
    
    for i in range(1, Nr-1):
    
        h = r[i+1] - r[i]  # assume approximately uniform
        
        # First derivative (central)
        dU[i] = (U[i+1] - U[i-1]) / (2*h)
        
        # Second derivative (central)
        d2U[i] = (U[i+1] - 2*U[i] + U[i-1]) / (h**2)
        
    # Right boundary (backward difference)
    
    h = r[-1] - r[-2]
    
    dU[-1] = (U[-1] - U[-2]) / h
    d2U[-1] = (U[-1] - 2*U[-2] + U[-3]) / (h**2)

    for i in range(0, Nr):
    
        D[i] = D0 * np.exp(-U[i]/(kB*T))
        dD[i] = -D0/(kB*T) * np.exp(-U[i]/(kB*T)) * dU[i]
        d2D[i] = (-D0/(kB*T)) * ( np.exp(-U[i]/(kB*T)) * d2U[i] + (-dU[i]/(kB*T)) * np.exp(-U[i]/(kB*T)) * dU[i] )

        mu[i] = D[i]/(kB*T)
        dmu[i] = dD[i]/(kB*T)

else:

    print("ERROR: Unknown potential definition flag: potential_def.")
    sys.exit(1)
    
##############################################################################
# Solving linear system
##############################################################################

def linear_system(lam):
    
    # Initialization
    
    A[:,:] = 0
    
    # Left boundary (forward difference)
    
    h = r[1] - r[0]
    
    A[0,0] = (2*r[0]**2*D[0]/(D0*h**2))
    A[0,1] = (r[0]**2*d2D[0]/D0) - (2*r[0]**2*D[0]/(D0*h**2)) + (2*r[0]*dD[0]/D0) + (2*r[0]*mu[0]*dU[0]/D0) + (r[0]**2*dmu[0]*dU[0]/D0) + (r[0]**2*mu[0]*d2U[0]/D0) + (lam**2*r[0]**2 - n*(n+1))
    
    # Local spacing approximation
    
    for j in range(1,Nr-1):
    
        h = r[j+1] - r[j]  # assume approximately uniform
        
        A[j,j-1] = (r[j]**2*D[j]/(D0*h**2)) - (r[j]**2*dD[j]/(D0*h)) - (r[j]*D[j]/(D0*h)) - (mu[j]*r[j]**2*dU[j]/(2*D0*h))
        A[j,j] = (r[j]**2*d2D[j]/D0) - (2*r[j]**2*D[j]/(D0*h**2)) + (2*r[j]*dD[j]/D0) + (2*r[j]*mu[j]*dU[j]/D0) + (r[j]**2*dmu[j]*dU[j]/D0) + (r[j]**2*mu[j]*d2U[j]/D0) + (lam**2*r[j]**2 - n*(n+1))
        A[j,j+1] = (r[j]**2*D[j]/(D0*h**2)) + (r[j]**2*dD[j]/(D0*h)) + (r[j]*D[j]/(D0*h)) + (mu[j]*r[j]**2*dU[j]/(2*D0*h))
    
    # Right boundary (backward difference)
    
    h = r[-2] - r[-1]
    
    A[Nr-1,Nr-2] = (2*r[-1]**2*D[-1]/(D0*h**2))
    A[Nr-1,Nr-1] = (r[-1]**2*d2D[-1]/D0) - (2*r[-1]**2*D[-1]/(D0*h**2)) + (2*r[-1]*dD[-1]/D0) + (2*r[-1]*mu[-1]*dU[-1]/D0) + (r[-1]**2*dmu[-1]*dU[-1]/D0) + (r[-1]**2*mu[-1]*d2U[-1]/D0) + (lam**2*r[-1]**2 - n*(n+1))
    
    # Solving linear system

    B[:] = 1e-60
    
    Rsi[:] = np.linalg.solve(A,B)
    
    Bsoli[:] = A[:,:]@Rsi[:]
    
    return Rsi,Bsoli


##############################################################################
# Finding the roots (eigenvalues) and tau's
##############################################################################

count = 0

for i in range(0,Nx):

    lam = x[i]

    Rsi,Bsoli = linear_system(lam)
    
    # Calculate the derivatives
    
    derivi = (Rsi[1] - Rsi[0])/(dr)
    derivf = (Rsi[-2] - Rsi[-1])/(dr)
    
    hi[i] = derivi
    hf[i] = derivf
    
    # Checking B.C.s (eigenvalues)

    if abs(hi[i]) + abs(hi[i-1]) != abs(hi[i]+hi[i-1]) and abs(hf[i]) + abs(hf[i-1]) != abs(hf[i]+hf[i-1]):

        count += 1
        
        if count>=1:
            
            idx.append(i)
            
            # Eigenvalue
            
            root.append( (x[i]+x[i-1])/2 )
            tau.append( 1/(root[-1]**2 * D0) )
            
            # Eigenfunction
            
            if count ==1:
                Rs = [[Rsi[j]] for j in range(len(Rsi))]
            else:
                Rs = [row + [Rsi[j]] for j, row in enumerate(Rs)]

            # Bsol vector
            
            if count ==1:
                Bsol = [[Bsoli[j]] for j in range(len(Bsoli))]
            else:
                Bsol = [row + [Bsoli[j]] for j, row in enumerate(Bsol)]

            # Printing the tau's

            print("Eigenvalue [√(1/m)] = "+"{:.10E}".format(root[-1]),"|","Tau [s] = "+"{:.10E}".format(tau[-1]))
                
            # Cut-off for the number of eigenmodes
            
            if count >= Km:
                
                break

##############################################################################
# Normalizing the radial solution
##############################################################################

lidx = len(idx)
Rs = np.array(Rs)
Bsol = np.array(Bsol)
r = np.append(r, r[-1]+dr)

def norm(Rs):
    
    I = 0
    
    for i in range(0,Nr):
        
        dr = r[i+1] - r[i]
        
        I += r[i]**2 * Rs[i]**2 * dr
    
    return I

# Normalization

for j in range(0,lidx):

    Rs[:,j] = Rs[:,j] / np.sqrt(norm(Rs[:,j]))
    
    
##############################################################################
# Calculating P(tau)
##############################################################################

Ptau = np.zeros([lidx],float)

# Calculating p(r0)

Nint = 0

for j in range(0,Nr):

    dr = r[j+1] - r[j]
    
    Nint += r[j]**2 * np.exp(-U[j]/(kB*T)) * dr

prefac = 1/(4*pi*Nint)

# Semi-normalized spherical harmonics

harmonics_norm = (16*pi)/5

# Calculating the constants

for i in range(0,lidx):
    
    Ptau[i] = 0
    
    for j in range(0,Nr): # r0 integration
    
        Arg = 0
        
        Ak = Rs[j,i]
    
        for k in range(0,Nr): # r integration
        
            dr = r[k+1] - r[k]
    
            Arg += Ak * Rs[k,i] / r[k] * dr
        
        dr = r[j+1] - r[j]
                
        Ptau[i] += alpha * prefac * harmonics_norm * np.exp(-U[j]/(kB*T)) * Arg / r[j] * dr
        
##############################################################################
# Calculating total G(t)
##############################################################################

for i in range(0,Nt):
    
    for j in range(0,lidx):
    
        Gt[i] += np.exp(-t[i]/tau[j]) * Ptau[j]
    
##############################################################################
# Calculating T12(0)
##############################################################################

T12i = 0

for j in range(0,lidx):
    
    T12i += 10 * Ptau[j] * tau[j]

print("1/T₁₂(0) [1/s] = ",T12i)

##############################################################################
# Exporting the data
##############################################################################

# Radial Solution (Eigenfunctions)

f = open("01_01_R_sol.dat","w")
sys.stdout = f

header = "r [m]\t" + "\t".join([f"R_{d+1}(r) [∛(1/m)] \t" for d in range(lidx)])
print(header)
for i in range(0,Nr):
    print("{:.6E}".format(r[i]), Rs[i,:])
                        
f.close()
sys.stdout = original_stdout

# Tau and P(tau)

f = open("01_02_Ptau.dat","w")
sys.stdout = f

print("tau [s]   P(tau) [1/s²]")
for i in range(0,lidx):
    print("{:.6E}".format(tau[i]),"{0:0.6f}".format(Ptau[i]))
                        
f.close()
sys.stdout = original_stdout

# NMR dipole-dipole autocorrelation function

f = open("01_03_Gt.dat","w")
sys.stdout = f

print("t [s]   G(t) [1/s²]")
for i in range(0,Nt):
    print("{:.6E}".format(t[i]),"{:.6E}".format(Gt[i]))
                        
f.close()
sys.stdout = original_stdout

# Relaxation rates at zero frequency limit

f = open("01_04_T12i.dat","w")
sys.stdout = f

print("1/T₁₂(0) [1/s]")
print(T12i)
                        
f.close()
sys.stdout = original_stdout
        
##############################################################################
# Plotting the results
##############################################################################

plt.figure(0)
for i in range(0,lidx):    
    plt.plot(r[:-1], Rs[:,i], color=plt.cm.viridis(i/lidx), linestyle='-', linewidth=1.0, label='lambda = '+str("{0:0.6f}".format(root[i])))
    plt.legend()
    plt.axhline(y=0, color='b', linestyle='--')
    plt.title('R(r) vs r')
    plt.ylabel('R(r) [∛(1/m)]')
    plt.xlabel('r [m]')
    plt.savefig("01_01_R_sol.png", dpi=500)
    
plt.figure(1)
for i in range(0,lidx):
    plt.vlines(tau[i], 0, Ptau[i], colors=plt.cm.viridis(i/lidx), linestyles='solid')
plt.xscale('log')
plt.title('P(τ) vs τ')
plt.ylabel('P(τ) [1/s²]')
plt.xlabel('τ [s]')
plt.savefig("01_02_Ptau.png", dpi=500)

plt.figure(2)
plt.plot(t[:], Gt[:], color='b', linestyle='-', linewidth=1.0)
plt.title('G(t) vs t')
plt.ylabel('G(t) [1/s²]')
plt.xlabel('t [s]')
plt.savefig("01_03_Gt.png", dpi=500)

plt.figure(3)
plt.plot(t[:], Gt[:], color='b', linestyle='-', linewidth=1.0)
plt.yscale('log')
plt.title('ln G(t) vs t')
plt.ylabel('G(t) [1/s²]')
plt.xlabel('t [s]')
plt.savefig("01_03_lnGt.png", dpi=500)
