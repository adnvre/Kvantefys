#!/usr/bin/env python
"""
Created on Sun 27 Oct 2013.
Modified to solve the home-exam in FYS2140 March 2015.


Uses the Crank Nicolson scheme to solve the time dependent Schrodinger equation
for a potential barrier.

Animation is done using the matplotlib.pyplot library.

@author Are Raklev
@author Benedicte Emilie Braekken
"""

# Tools for sparse matrices
import scipy.sparse as sparse
import scipy.sparse.linalg

# Numerical tools
from numpy import *

# Plotting library
from matplotlib.pyplot import *

"""Physical constants"""
_E0e =   0.511       # Rest energy for an electron [MeV]
_hbarc = 0.1973      # [MeV pm]
_c = 3.0e2           # Spees of light [pm / as]

"""Parameters of initial state"""
a  =  1.00 # [pm]
l  =  2.0 # [1 / pm]

# Define functions
def Psi0( x , a, l ):
    '''
    Initial state for a travelling gaussian wave packet.
    '''
    x0 = -5.00 # [pm]
    
    A = ( 1. / ( 2 * pi * a**2 ) )**0.25
    K1 = exp( - ( x - x0 )**2 / ( 4. * a**2 ) )
    K2 = exp( 1j * l * x )
                  
    return A * K1 * K2#5

def potentialBarrier( x, height=1, width=0.0025 ):
    """
    Gives the potential for a potential well of depth depth and width width.

    @param height Gives the height of the potential well. Given as the magnitude
    (positive integer / double / float).
    @param width Gives the width of the potential well. Given as positive
    integer definig the fraction of the x spectrum to contain the well. For
    example, 1 will mean that the well covers the whole spectrum and 0.5 that
    it covers half of it.
    """
    # Declare new empty array with same length as x
    potential = zeros( len( x ) )

    potential[ 0.5*len(potential) : (0.5+width)*len(potential) ] = height

    return potential

if __name__ == '__main__':
    nx = 5001 # Number of points in x direction
    dx = 0.020 # Distance between x points [pm]

    # Use zero as center, same amount of points each side
    x1 = - 0.5 * nx * dx
    x2 = 0.5 * nx * dx
    x = linspace( x1, x2, nx )

    # Time parameters
    T = 0.100 # How long to run simulation [as]
    dt = 5e-4 # The time step [as]
    t = 0
    time_steps = int( T / dt ) # Number of time steps

    # Constants - save time by calculating outside of loop
    k1 = ( 1j * _hbarc * _c) / (2. * _E0e )
    k2 = - ( 1j * _c ) / _hbarc

    # Create the initial state Psi
    Psi = Psi0(x,a,l)

    # Create the matrix containing central differences. It it used to
    # approximate the second derivative.
    data = ones((3, nx))
    data[1] = -2*data[1]
    diags = [-1,0,1]
    D2 = k1 / dx**2 * sparse.spdiags(data,diags,nx,nx)

    # Identity Matrix
    I = sparse.identity(nx)

    # Create the diagonal matrix containing the potential.
    V_data = potentialBarrier(x)
    V_diags = [0]
    V = k2 * sparse.spdiags(V_data, V_diags, nx, nx)

    # Put mmatplotlib in interactive mode for animation
    ion()

    # Setup the figure before starting animation
    fig = figure() # Create window
    ax = fig.add_subplot(111) # Add axes
    line, = ax.plot( x, abs(Psi)**2, label='$|\Psi(x,t)|^2$' ) # Fetch the line object

    # Also draw a green line illustrating the potential
    ax.plot( x, V_data, label='$V(x)$' )

    # Add other properties to the plot to make it more elegant
    fig.suptitle("Solution of Schrodinger's equation with potential barrier") # Title of plot
    ax.grid('on') # Square grid lines in plot
    ax.set_xlabel('$x$ [pm]') # X label of axes
    ax.set_ylabel('$|\Psi(x, t)|^2$ [1/pm] and $V(x)$ [MeV]') # Y label of axes
    ax.set_xlim([-30, 30])  # Sets x-axis range
    ax.set_ylim([0, 1.1])   # Sets y-axis range
    ax.legend(loc='best')   # Adds labels of the lines to the window
    draw() # Draws first window

    # Time loop
    while t < T:
        """
        For each iteration: Solve the system of linear equations:
        (I - k/2*D2) u_new = (I + k/2*D2)*u_old
        """
        # Set the elements of the equation
    	A = (I - dt/2. * (D2 + V))
    	b = (I + dt/2. * (D2 + V)) * Psi

        # Calculate the new Psi
    	Psi = sparse.linalg.spsolve(A,b)

        # Update time
    	t += dt

    	# Plot this new state
    	line.set_ydata( abs(Psi)**2 ) # Update the y values of the Psi line
        draw() # Update the plot


    # Integral to find transmission probability for a given energy
    E = _hbarc**2*(l**2 + 1./(4.*a**2))/(2*_E0e)  # Calculate energy expectation value [MeV]
    Psi_pluss = Psi[(nx-1)/2:]                    # Slice away list for x < 0
    I = trapz(abs(Psi_pluss)**2,None,dx)          # Integrate with trapezoidal method
    print E, I

  
    # Turn off interactive mode
    ioff()

    # Add show so that windows do not automatically close
    show()