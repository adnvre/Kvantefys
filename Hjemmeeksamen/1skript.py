#!/usr/bin/env python
"""
Created on Sun 27 Oct 2013.
Modified to solve the home-exam in FYS2140 March 2015.

@author Benedicte Emilie Braekken
@author Are Raklev
"""

# Numerical tools
from numpy import *

# Plotting library
from matplotlib.pyplot import *

# Function definitions
def Psi0( x ):
    '''
    Initial state for a gaussian wave packet.
    '''
    x1 = -5.00   # [pm] Start here
    a  =  1.00   # [pm] Width of packet

    A = ( 1. / ( 2 * pi * a**2 ) )**0.25
    K1 = exp( - ( x - x1 )**2 / ( 4. * a**2 ) )

    return A * K1

# Main programme
if __name__ == '__main__':
    
    # Define some numbers for x-axis
    nx = 4001    # Number of points in x direction
    dx = 0.005 # Distance between x points [pm]

    # Interval [a,b], same amount of points each side
    a = - 0.5 * nx * dx
    b = 0.5 * nx * dx
    x = linspace( a, b, nx )   # x is now an array containing all x-values used
    
    # Plot initial state
    figure()        # Create a window for figure
    Psi = Psi0(x)   # Create the initial state as an array Psi from the array x
    plot( x, abs(Psi)**2, label='$|\Psi(x,t)|^2$' ) # Plot
    xlabel('$x$ [pm]')                # Label for x-axis
    ylabel('$|\Psi(x, t)|^2$ [1/pm]') # Label for y-axis
    legend(loc='best')                # Adds labels of the lines to the window
    savefig('psisq_init.eps')         # Save as .eps figure
    axis([-10, 10, 0, 0.7])           # Set axis range


    # Turn off interactive mode
    ioff()

    # Add show so that windows do not automatically close
    show()