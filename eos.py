########################################################################
# Team Rabbits: Austin Czyzewski | Megan Davis | Molly Janasik | David Pham
# AST304, Fall 2018
# Michigan State University
#   
#
# The purpose of this module is to define a density and pressure relationship so that we can use either form of data
# in structure.py 

########################################################################

# import/define constants here
import astro_const as ac

def pressure(rho, mue):
    """
    Arguments
        rho
            mass density (kg/m**3)
        mue
            baryon/electron ratio
    
    Returns
        electron degeneracy pressure (Pascal)
    """
    
    #p = ((ac.pi**2*ac.hbar**2)/(5*ac.m_e*ac.m_p**(5/3)))*(3/ac.pi)**(2/3)*(rho/mue)**(5/3)
    
    # from equation 1 from GCP2 project file
    p = (1/5)*(3/(8*ac.pi))**(2/3)*(ac.h)**2/ac.m_e*(rho/(mue*ac.m_u))**(5/3) 
    
    return p 

def density(p, mue):
    """
    Arguments
        p
            electron degeneracy pressure (Pascal)
        mue
            baryon/electron ratio
        
    Returns
        mass density (kg/m**3)
    """
    
    #rho = (mue)*(((5*ac.m_e*ac.m_p**(5/3)*p)/(ac.pi**2*ac.hbar**2))*(ac.pi/3)**(2/3))**(3/5)
    
    # from equation 1 from GCP2 project file
    rho = ((5*p)*((8*ac.pi)/3)**(2/3)*((ac.m_e)/(ac.h**2)))**(3/5)*(mue*ac.m_u) 
    
    return rho