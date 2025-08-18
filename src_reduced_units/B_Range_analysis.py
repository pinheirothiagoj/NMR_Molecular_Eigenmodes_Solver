##############################################################################
# RICE UNIVERSITY
# Department of Chemical and Biomolecular Engineering
# Author: Thiago Jose Pinheiro dos Santos, Ph.D.
##############################################################################

# This code computes the multiexponential time decays (τ's), their amplitudes
# P(τ)'s, and 1/T1(0) = 1/T2(0) for different ranges of ri/rf, with or without
# a potential of interaction. This code reads input data from 'input.txt', and
# it also modifies this input file to update ri and rf values accordingly.
# This code (B) is intended to be used in conjunction with code (A).
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
import os
import sys

###############################################################################
# Setting up variables
###############################################################################

# Spherical thickness discretization

dri = 0.025
ri = np.arange(0.40,2.00,dri)
rf = 2.00

# Number of eigenmodes

Km = 5

###############################################################################
# Defining vectors and matrices
###############################################################################

Nlen = len(ri)

tau = np.zeros([Nlen,Km],float)
Ptau = np.zeros([Nlen,Km],float)
T12i = np.zeros([Nlen],float)

original_stdout = sys.stdout
np.set_printoptions(threshold=sys.maxsize)

##############################################################################
# Loop over different ri/rf
##############################################################################

for i in range(0,Nlen):

    # Rewriting the input file
    
    filename = "input.txt"
    replacements = {'ri': ri[i], 'rf': rf, 'Km': Km}
    replaced = set()
    
    with open("input.txt", 'r') as file:
        lines = []
        for line in file:
            for key, value in replacements.items():
                if key + ' = ' in line and key not in replaced:
                    line = f"{key} = {value}\n"
                    replaced.add(key)
            lines.append(line)
    
    with open("input.txt", 'w') as file:
        file.writelines(lines)
        
    # Calling the numerical solver

    bashCommand = "python ../../src_reduced_units/A_Numeric_radial_FP_reducedunits.py"
    print("Solving for ri/rf = ", ri[i]/rf)
    os.system(bashCommand)

    f = open("01_02_Ptau.dat","r")
    content = f.readlines()[1:]
    
    l = len(content)
    
    for d in range(0,l):
        
        T = content[d]

        tau[i,d] = T.split()[0]   # tau
        Ptau[i,d] = T.split()[1]  # P(tau)
        
    f.close()

    f = open("01_04_T12i.dat","r")
    content = f.readlines()[1:]
    
    T = content[0]
    T12i[i] = T.split()[0]
        
    f.close()

###############################################################################
# Exporting the data
###############################################################################

# Tau and P(tau) in different radial shells

f = open("02_01_profile_Ptau.dat","w")
sys.stdout = f

header = "ri/rf\t" + "\t".join([f"Mode_{d+1}" for d in range(Km)])
print(header)
for i in range(0,Nlen):
    print("{0:0.6f}".format(ri[i]/rf), *["{0:0.6f}".format(val) for val in Ptau[i,:]], sep="\t")
                        
f.close()
sys.stdout = original_stdout

# Relaxation rates at zero frequency in different radial shells

f = open("02_02_profile_T12i.dat","w")
sys.stdout = f

print("ri/rf\t1/T₁₂(0)")
for i in range(0,Nlen):
    print("{0:0.6f}".format(ri[i]/rf), "{0:0.6f}".format(T12i[i]), sep="\t")
                        
f.close()
sys.stdout = original_stdout


###############################################################################
# Plotting the results
###############################################################################

plt.figure(0)   
for d in range(0,Km): 
    plt.plot(ri[:]/rf, Ptau[:,d], color=plt.cm.viridis(d/Km), linestyle='-', linewidth=1.0, label='Mode = '+str(d+1))
    plt.legend()
    plt.axhline(y=0, color='b', linestyle='--')
    plt.title('P(τ) vs ri/rf')
    plt.ylabel('P(τ)')
    plt.xlabel('ri/rf')
    plt.savefig("02_01_profile_Ptau.png", dpi=500)

plt.figure(1)    
for d in range(0,Km):    
    plt.plot(ri[:]/rf, tau[:,d], color=plt.cm.viridis(d/Km), linestyle='-', linewidth=1.0, label='Mode = '+str(d+1))
    plt.legend()
    plt.axhline(y=0, color='b', linestyle='--')
    plt.title('τ vs ri/rf')
    plt.ylabel('τ')
    plt.xlabel('ri/rf')
    plt.savefig("02_01_profile_tau.png", dpi=500)

plt.figure(2)    
plt.plot(ri[:]/rf, T12i[:], color='b', linestyle='-', linewidth=1.0)
plt.title('1/T₁₂(0) vs ri/rf')
plt.ylabel('1/T₁₂(0)')
plt.xlabel('ri/rf')
plt.savefig("02_02_profile_T12i.png", dpi=500)
