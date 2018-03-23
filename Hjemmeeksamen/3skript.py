from numpy import *
from matplotlib.pyplot import *


"""Physical constants"""
_E0e =   0.511       # Rest energy for an electron [MeV]
_hbarc = 0.1973      # [MeV pm]


# Function definition
def T( E ):
    '''
    Expression for transmission probability
    '''
    V0 = 34.   # Height of potential [MeV]
    L = 17.  # Length of potential box [pm]
    
    T = 1. / (1. + V0**2/(4.*(V0-E)*E)*sinh(sqrt(2.*_E0e*(V0-E))/_hbarc*L)**2)

    return T

# Main programme
if __name__ == '__main__':
    
    # Define some numbers for axis
    nE = 400    # Number of points
    dE = 0.0025  # Distance between points [MeV]

    # Interval [a,b]
    a = 0
    b = nE * dE
    E = linspace( a, b, nE )   # An array containing all values of energy used
    
    # Plotting
    figure()           # Create a window for figure
    plot( E, T(E), label='$T(E)$ analytisk' ) # Plot T(E) function
    Evalues = [0.009,0.019,0.048,0.162,0.352,0.619,0.962,1.381]
    Tvalues = [0.005,0.022,0.057,0.181,0.348,0.520,0.671,0.789]
    plot(Evalues, Tvalues, 'ro', label='$T(E)$ numerisk' ) # Plot data points
    xlabel('$E$ [MeV]')               # Label for x-axis
    ylabel('$T(E)$')                  # Label for y-axis
    legend(loc='best')                # Adds labels of the lines to the window
    savefig('TE1.eps')                # Save as .eps figure
    axis([0, 1.5, 0, 1.0])              # Set axis range

    # Plotting
    figure()           # Create a window for figure
    plot( E, T(E), label='$T(E)$ analytisk' ) # Plot T(E) function
    plot(Evalues, Tvalues, 'ro', label='$T(E)$ numerisk') # Plot data points
    Evalues = [0.009,0.018,0.044,0.083,0.130,0.179]
    Tvalues = [0.028,0.090,0.177,0.257,0.315,0.343]
    plot(Evalues, Tvalues, 'rx',label='$T(E)$ rel. korr.') # Plot data points
    xlabel('$E$ [MeV]')               # Label for x-axis
    ylabel('$T(E)$')                  # Label for y-axis
    legend(loc='best')                # Adds labels of the lines to the window
    axis([0, 0.3, 0, .5])             # Set axis range
    savefig('TE2.eps')                # Save as .eps figure


    # Turn off interactive mode
    ioff()

    # Add show so that windows do not automatically close
    show()