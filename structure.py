########################################################################
# Team Rabbits: Austin Czyzewski | Megan Davis | Molly Janasik | David Pham
# AST304, Fall 2018
# Michigan State University
#   
# This module will set the central boundary condition and construct the integration loop.
########################################################################

import numpy as np
#from numpy.linalg import norm
from eos import density
from ode import rk4
import astro_const as sc


def stellar_derivatives(m,z,mue):
    """
    RHS of Lagrangian differential equations for radius and pressure
    
    Arguments
        m
            current value of the mass (units = kg)
        z (array)
            current values of (radius, pressure) (units = m, Pa)
        
    Returns
        dzdm (array)
            Lagrangian derivatives dr/dm, dP/dm (units =m/kg,Pa/kg)
    """
    rho = density(z[1],mue)
    drdm = 1/(4*np.pi*z[0]**2*rho)  #Lagrangian Derivative of r w.r.t. m
    dpdm = (-sc.G*m) /(4 * z[0]**4 * np.pi) #Lagrangian Derivative of p w.r.t. m
    
    dzdm = np.zeros_like(z)
    dzdm[0] = drdm
    dzdm[1] = dpdm
    

    return dzdm

def central_values(Pc,delta_m,mue):
    """
    Constructs the boundary conditions at the edge of a small, constant density 
    core of mass delta_m with central pressure P_c
    
    Arguments
        Pc
            central pressure (units = Pa)
        delta_m
            core mass (units = kg)
        mue
            nucleon/electron ratio
    
    Returns
        z = array([ r, p ])
            central values of radius and pressure (units = m, Pa)
    """
    z = np.zeros(2)
    
    r_delta_m = ((3*delta_m) / (4*sc.pi*density(Pc, mue)))**(1/3) #Radius depends on delta_m and density
    p_delta_m = Pc #Pressure is equal to Pc
    
    # compute initial values of z = [ r, p ]
    z = [r_delta_m , p_delta_m]
    
    return z

def lengthscales(m,z,mue):
    """
    Computes the radial length scale H_r and the pressure length H_P
    
    Arguments
        m
            current mass coordinate (units = kg)
        z (array)
           [ r, p ] (units = m, Pa)
        mue
            mean electron weight
    
    Returns
        z/|dzdm| (units = kg,kg)
    """

    dzdm = stellar_derivatives(m,z,mue)
    Hr = z[0] / np.abs(dzdm[0])
    Hp = z[1] / np.abs(dzdm[1])
    zdzdm = [Hr,Hp]
    
    return zdzdm

    ##########
    # Old code that does the same thing as what is in the above lines
    ##########
    # z/abs(stellar_derivatives(m,z,mue))
    
    #return  

    #Since we are simply dividing z by the Langrangian derivatives we can
    #do this in one line given the inputs
    
def integrate(Pc,delta_m,eta,xi,mue,max_steps=10000):
    """
    Integrates the scaled stellar structure equations

    Arguments
        Pc
            central pressure (units = Pa)
        delta_m
            initial offset from center (units = kg)
        eta
            The integration stops when P < eta * Pc
        xi
            The stepsize is set to be xi*min(p/|dp/dm|, r/|dr/dm|)
        mue
            mean electron mass
        max_steps
            solver will quit and throw error if this more than max_steps are 
            required (default is 10000)
                        
    Returns
        m_step, r_step, p_step
            arrays containing mass coordinates, radii and pressures during 
            integration (kg,m,Pa)
    """
    
    # making arrays of zeros
    m_step = np.zeros(max_steps)
    r_step = np.zeros(max_steps)
    p_step = np.zeros(max_steps)
    
    # initial conditions
    
    z = central_values(Pc,delta_m,mue)
    m_step[0] = delta_m
    r_step[0] = z[0]
    p_step[0] = z[1]
    
    zdzdm = lengthscales(delta_m,z,mue)
    
    #initial stepsize set
    h = xi * min(zdzdm[0],zdzdm[1])
    
    m = delta_m
#    advance_one_step = integration_methods[rk4]
    
    Nsteps = 0
    for step in range(1,max_steps):
        
    
        # check for completion
        if (z[1] < eta*Pc):
            break
        
        h = xi * min(lengthscales(m,z,mue))
        # store the step
        
        z = rk4(stellar_derivatives,m,z,h,args=mue)
        
        m_step[step] = m
        r_step[step] = z[0]
        p_step[step] = z[1]
        
        # set the stepsize
        #zdzdm = lengthscales(m,z,mue)
       
        
        # take a step
        #z = advance_one_step(stellar_derivatives,z[0],z,h,rk4)
        
        #increment the counter
        Nsteps += 1
        m += h
    # if the loop runs to max_steps, then produce an error
    else:
        raise Exception('too many iterations')
        
    return m_step[0:Nsteps],r_step[0:Nsteps],p_step[0:Nsteps]


def pressure_guess(m,mue):
    """
    Computes a guess for the central pressure based on virial theorem and mass-
    radius relation. 
    
    Arguments
        m
            mass of white dwarf (units = kg)
        mue
            mean electron mass (units = kg)
    
    Returns
        P
            guess for pressure (units = Pa)
    """
    
    Pguess = (sc.G**5/sc.Ke**4)*(m*mue**2)**(10/3) #use imported constants, using equation 16
    
    return Pguess


# old stuff we do not want to delete
#delta_m = sc.Msun * 1e-10
#print(np.shape(integrate(pressure_guess(1e30,2),1e30,10e-10,0.1,sc.m_e)))
#plt.plot(integrate(pressure_guess(sc.Msun,2),delta_m,10e-10,0.01,2))
#plt.yscale('log')
#plt.show()
