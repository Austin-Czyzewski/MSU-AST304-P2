########################################################################
# Team Rabbits: Austin Czyzewski | Megan Davis | Molly Janasik | David Pham
# AST304, Fall 2018
# Michigan State University
#   
# Useful astrophysical constants
# Solar values were taken from the Particle Data Group,
# http://pdg.lbl.gov/2018/reviews/rpp2018-rev-astrophysical-constants.pdf
# The remaining are pulled in from scipy.constants
########################################################################

import scipy.constants as _sc
import numpy as np

# solar mass, radius, luminosity
Msun = 1.98848e30 # kg
Rsun = 6.957e8 # m
Lsun = 3.828e26 # W

# physical constants from scipy, all in MKS units
G = _sc.G
h = _sc.h
hbar = _sc.hbar
m_e = _sc.m_e
m_p = _sc.m_p
m_n = _sc.m_n
m_u = _sc.m_u
c = _sc.c
k = _sc.k
pc = _sc.parsec
au = _sc.astronomical_unit
year = _sc.year
sigmaSB = _sc.sigma
pi = _sc.pi
Ke = (1/5) * (3/(8*pi))**(2/3) * h**2 / m_e * (1/(m_u*2))**(5/3) #defined constant for structure

def mean_molecular_weight(n,m):
    '''
    n (array-like)
        Array of the values of percentage of composition or total value of compostion
    m (Array-like) 
        In the same order, arrays of the mass of each molecule
    
    returns
        mean molecular weight (units = kg)
    '''
    weight_total = 0
    number_total = 0
    for i in range(len(n)):
        weight_total += n[i]*m[i]
        number_total += n[i]
    return weight_total/number_total

Mass_H = m_p+m_e #in kg
Mass_He = 6.6465e-27 #kg
Mass_N = 1.6749e-27 #kg
Num_H = 0.706
Num_He = 0.275
Num_N = 0.019
Number_Densities = np.array([Num_H,Num_He,Num_N])
Molecular_Weights = np.array([Mass_H,Mass_He,Mass_N])
mue_interstellar = mean_molecular_weight(Number_Densities,Molecular_Weights)

print("The Mean Molecular Weight for interstellar gas is {:.5}".format(mue_interstellar))

if __name__ == "__main__":
    
    constants = [
        ("solar mass",Msun,"kg"),
        ("solar radius",Rsun,"m"),
        ("solar luminosity",Lsun,"W"),
        ("gravitational constant",G,"m**3 s**-2 kg**-1"),
        ("Planck constant",h,"J s"),
        ("Planck constant, reduced",hbar,"J s"),
        ("electron mass",m_e,"kg"),
        ("proton mass",m_p,"kg"),
        ("neutron mass",m_n,"kg"),
        ("atomic mass unit",m_u,"kg"),
        ("speed of light",c,"m s**-1"),
        ("Boltzmann constant",k,"J K**-1"),
        ("parsec",pc,"m"),
        ("astronomical unit",au,"m"),
        ("year",year,"s"),
        ("Stefan-Boltzmann constant",sigmaSB,"W m**-2 K**-4")
    ]
    
    
    for const in constants:
        print('{0[0]:28} = {0[1]:11.4e} {0[2]}'.format(const))