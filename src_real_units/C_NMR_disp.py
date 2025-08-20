##############################################################################
# RICE UNIVERSITY
# Department of Chemical and Biomolecular Engineering
# Author: Thiago Jose Pinheiro dos Santos, Ph.D.
##############################################################################

# This code computes the dipole-dipole NMR frequency dispersion for 1/T1 and
# 1/T2, for a single case of dipole pairs with either like or unlike spins.
# This code reads input data from 'input.txt'. This code (C) is intended to be
# used in conjunction with code (A)
#
# When using these codes or any of their parts, please cite:
#
# [1] T. J. Pinheiro dos Santos, B. Orcan-Ekmecki, W. G. Chapman, P. M. Singer,
# and D. N. Asthagiri, “Extended molecular eigenmodes treatment of dipole-
# dipole NMR relaxation in real fluids”.
#
# [2] T. J. Pinheiro dos Santos, B. Orcan-Ekmecki, W. G. Chapman, P. M. Singer,
# and D. N. Asthagiri, “Theory and modeling of molecular modes in the NMR
# relaxation of fluids,” J. Chem. Phys., vol. 064108, p. 160, 2024.


import numpy as np
import matplotlib.pyplot as plt
import os
import sys

##############################################################################
# Setting up variables
##############################################################################

# Module for like spins ("ls") or unlike spins ("us")

mod_spin = "ls"

# Execute code A beforehand ("yes" or "no")

exeA = "yes"

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

# Frequency discretization

Nw = int(1e5)                                                   # units = Hz
w0 = np.logspace(np.log10(1e3),np.log10(1e10),num=Nw,base=10)   # units = Hz

##############################################################################
# Defining vectors and matrices
##############################################################################

Nlen = len(w0)

tau = np.zeros([Km],float)
Ptau = np.zeros([Km],float)
T1i = np.zeros([Nlen],float)
T2i = np.zeros([Nlen],float)

original_stdout = sys.stdout
np.set_printoptions(threshold=sys.maxsize)

##############################################################################
# Calling the numerical solver
##############################################################################

if exeA == "yes":
        
    bashCommand = "python ../../src_real_units/A_Numeric_radial_FP_realunits.py"
    print("Solving for ri/rf = ", ri/rf)
    os.system(bashCommand)

elif exeA == "no":

    print("WARNING: Skipping the execution of code A and reading τ and P(τ) from the file.")

else:

    print("ERROR: Unknown module A executor flag: exeA.")
    sys.exit(1)

# Reading the tau's and P(tau)'s

f = open("01_02_Ptau.dat","r")
content = f.readlines()[1:] 

l = len(content)

for d in range(0,Km):
    
    Tl = content[d]

    tau[d] = Tl.split()[0]   # tau
    Ptau[d] = Tl.split()[1]  # P(tau)
    
f.close()

##############################################################################
# Calculating NMR dispersion
##############################################################################

for i in range(0,Nlen):
    
    T1i[i] = 0
    T2i[i] = 0
    
    if mod_spin == "ls":
        
        # Case of like spins
    
        for k in range(0,Km):
            
            T1i[i] += 2 * Ptau[k] * ( tau[k]/(1+(w0[i]*tau[k])**2) + 4*tau[k]/(1+(2*w0[i]*tau[k])**2) )
            T2i[i] += 2 * Ptau[k] * ( (3/2)*tau[k]/(1+(0*tau[k])**2) + (5/2)*tau[k]/(1+(w0[i]*tau[k])**2) + tau[k]/(1+(2*w0[i]*tau[k])**2) )

    elif mod_spin == "us":
        
        ws = 658 * w0[i]   # units = Hz
        
        # Case of unlike spins
    
        for k in range(0,Km):
            
            T1i[i] += 2 * Ptau[k] * ( (1/3)*tau[k]/(1+((w0[i]-ws)*tau[k])**2) + tau[k]/(1+(w0[i]*tau[k])**2) + 2*tau[k]/(1+((w0[i]+ws)*tau[k])**2) )
            T2i[i] += 2 * Ptau[k] * ( (2/3)*tau[k]/(1+(0*tau[k])**2) + (1/6)*tau[k]/(1+((w0[i]-ws)*tau[k])**2) + (1/2)*tau[k]/(1+(w0[i]*tau[k])**2) + tau[k]/(1+(ws*tau[k])**2) + tau[k]/(1+((w0[i]+ws)*tau[k])**2) )
            
    else:

        print("ERROR: Unknown module for like- or unlike-spins flag: mod_spin.")
        
        sys.exit(1)


##############################################################################
# Exporting the data
##############################################################################

# NMR dispersion of the longitudinal relaxation rate

f = open("03_01_T1i_disp.dat","w")
sys.stdout = f

print("ω₀ [Hz]      1/T₁ [1/s]")
for i in range(0,Nlen):
    print("{0:0.6f}".format(w0[i]), T1i[i])
                        
f.close()
sys.stdout = original_stdout

# NMR dispersion of the transverse relaxation rate

f = open("03_01_T2i_disp.dat","w")
sys.stdout = f

print("ω₀ [Hz]      1/T₂ [1/s]")
for i in range(0,Nlen):
    print("{0:0.6f}".format(w0[i]), T2i[i])
                        
f.close()
sys.stdout = original_stdout

##############################################################################
# Plotting the results
##############################################################################

plt.figure(0)    
plt.plot(w0, T1i, color='b', linestyle='-', linewidth=1.0)
plt.xscale('log')
plt.yscale('log')
plt.title('1/T₁(ω₀) vs ω₀')
plt.ylabel('1/T₁(ω₀) [1/s]')
plt.xlabel('ω₀ [Hz]')
plt.savefig("03_01_T1i_disp.png", dpi=500)

plt.figure(1)    
plt.plot(w0, T2i, color='b', linestyle='-', linewidth=1.0)
plt.xscale('log')
plt.yscale('log')
plt.title('1/T₂(ω₀) vs ω₀')
plt.ylabel('1/T₂(ω₀) [1/s]')
plt.xlabel('ω₀ [Hz]')
plt.savefig("03_02_T2i_disp.png", dpi=500)
