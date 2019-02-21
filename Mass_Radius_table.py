########################################################################
# Team Rabbits: Austin Czyzewski | Megan Davis | Molly Janasik | David Pham
# AST 304, Fall 2018
# Michigan State University
########################################################################
import numpy as np
from eos import density, pressure
from structure import integrate, pressure_guess
import astro_const as sc
from scipy.optimize import brentq
#######################################################################
# Setting up initial values for inputs
#######################################################################
xi = .01
eta = 10**-10
mue = 2
delta_m = sc.Msun * 1e-10
actual_masses = np.arange(.1,1.1,.1) * sc.Msun

# creating the first pressure guess to find first masses, radii, densities
Pguess = pressure_guess(sc.Msun,mue)

#######################################################################
# Doing integration for Pressure and making arrays
#######################################################################

def mass_diff(Pguess, M_want, delta_m, eta, xi, mue):
    """
    Computes in 2.7 the difference between the wanted mass 
    and given central pressure guess
    
    Pguess: given central pressure guess
    M_want: wanted mass
    eta: gives when to stop integrating
    xi: part of the stepsize
    mue: mean electron mass
    
    """
    ms_pg, rs_pg, ps_pg = integrate(Pguess,delta_m,eta,xi,mue,max_steps=10000)
    difference = (ms_pg[-1] - M_want)
    
    return difference
    
#######################################################################
# Bisect function and creating the table!
#######################################################################
def get_zeros(mass):
    """
    Takes a low and high pressure, finds the actual central pressure
    Uses brentq to find zero using mass_diff as the function
    
    returns Pc (central pressure), ms (masses), rs (radii), ps (pressures)
    """
    # creating values for high and low pressure
    Low_pressure = 1e19
    High_pressure = 5e22
    
    #Get Central Pressure
    Pc = brentq(mass_diff,Low_pressure,High_pressure,args=(mass,delta_m,eta,xi,mue))
            
    #getting masses, radii, pressures
    ms, rs, ps = integrate(Pc, delta_m, eta, xi, mue)
    
    return Pc, ms, rs, ps
   
#######################################################################
# Let's create that table
#######################################################################
def create_table():
    """
    Creates tables of values and saves them as a .txt file
            
    returns: nothing explicitly, function creates table
    """
    outfile = 'table.txt'
    with open(outfile, 'w') as fout:
        fout.write(' | {0:7} | {1:7} | {2:7} | {3:7} | {4:7} | {5:7} | \n'.
                      format('M/M_sun ','R/R_sun','P_c (MKS)','P_c / GM^2R^-4','rho_c (MKS)','rho / [3M /4 *pi R^3] '))
        # '--:' means right-align
        fout.write(' | {0:7} | {1:7} | {2:7} | {3:7} | {4:7} | {5:7} | \n'.format('--:','--:','--:','--:','--:','--:'))
            
        for i in range(len(actual_masses)):
            mass = actual_masses[i]
                
            Pc, ms, rs, ps = get_zeros(mass)
                
            #getting actual central density and radius
            central_r = rs[-1]
            central_dens = density(ps[0],mue)
            
            fout.write(' | {0:7.1f} | {1:7e} | {2:7e} | {3:7e} | {4:7e} | {5:7f} | \n'.
                       format(mass/sc.Msun, central_r/sc.Rsun , Pc, Pc * central_r**4/(sc.G*mass**2), central_dens, 
                              central_dens / (3*mass/(4*np.pi*central_r**3))))

if __name__ == "__main__":
    print("Hey, yeah sorry, this might take a while to run")
    create_table()
         
        
    
    
    
    
    
    