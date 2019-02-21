########################################################################
# Team Rabbits: Austin Czyzewski | Megan Davis | Molly Janasik | David Pham
# AST 304, Fall 2018
# Michigan State University

#This script will create the Mass-Radius Relationship Plot and the Observation data
########################################################################

import matplotlib.pyplot as plt
import observations as o
import astro_const as sc
from structure import pressure_guess
import Mass_Radius_table as mrt
#######################################################################
# Setting up initial values for inputs
#######################################################################

radii=[]

#Collect radii

Pguess = pressure_guess(sc.Msun, mrt.mue)
for i in mrt.actual_masses:
    Pc, ms, rs, ps = mrt.get_zeros(i)
    radii.append(rs[-1]/(sc.Rsun))
    
#Get observation data from observation.py and Joyce.txt    
obs= o.MassRadiusObservations()


def create_fig(ms,rs):
    """
    Creates fig comparing calculated mass-radius relationship and observation data
    
    input:
        MUST BE NORMALIZED TO SUN
        ms: masses
        rs: radii
        
        
    returns: nothing explicitly, function creates plot
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #making it so that we can get the sun symbol and subscript
    params = {'mathtext.default': 'regular', 'axes.labelsize': 'x-large', 'xtick.labelsize':'x-large'}          
    plt.rcParams.update(params)
    ax.set_xlabel('$M/M_\u2609$')
    ax.set_ylabel('$R/R_\u2609$')
    ax.set_title('Mass-Radius Relationship')
    ax.scatter(ms, rs, marker='o', label='Calculated Data')
    ax.errorbar(obs.mass,(0.01)*obs.radius,\
            yerr=(0.01)*obs.radius_err,xerr=obs.mass_err,fmt='none', color= 'grey',markersize=4, alpha=0.75)
    ax.scatter(obs.mass,(0.01)*obs.radius, label='Joyce et al Data')
    ax.legend()
    ax.grid()
    plt.savefig('RvM_plot.jpg', pad_inches=2, overwrite=True)
    plt.show()
 
    
if __name__ == "__main__":
    print("hey, might take a second to run")
    create_fig(mrt.actual_masses/(sc.Msun), radii)










