# Numerical inverse Laplace transformation

# Integration
######################################
import scipy
from scipy import *
from scipy import integrate

# func: the function that you want to find its inverse
#
# t: the time step to calculate the inverse function at
#
# a: any arbitrary real point that func ...
# has no singularities to the right of it in the complex plane
#
# relError_tolerance: The tolerance for accuracy
#
# echo_out_process: A flag; set it to non-zero values to echo out the mid results.
def inverseLaplace(func, t, a, echo_out_process = 0, relError_tolerance = 1e-8):
    boundSteps = 5.0 / t # The initial segments of integration

    # Defining the bounds of integral
    lowerBound = 0
    upperBound = lowerBound + boundSteps

    # Defining the integrand based on argument function
    integF    = lambda u: real(func(complex(a+1j*u))) * cos(u*t)

    # Running the first integration outside the loop
    lastInegIter = 0
    integIter    = (2.0*exp(a*t) / pi) * integrate.quad(integF, lowerBound, upperBound)[0]

    segmentRuns = 1
    av_zero = 1e-15 # A very small number, to avoid devition by 0
    relError    = abs(abs(integIter) - abs(lastInegIter)) / abs(lastInegIter+av_zero)

    # Iteratively running the integration
    while relError > relError_tolerance:
        segmentRuns  = segmentRuns + 1
        lowerBound   = upperBound
        upperBound   = upperBound + boundSteps
        lastInegIter = integIter
        integIter    = integIter + (2.0*exp(a*t) / pi) * integrate.quad(integF, lowerBound, upperBound)[0]
        relError     = abs(abs(integIter) - abs(lastInegIter)) / (abs(lastInegIter)+av_zero)

    if (echo_out_process == 0):
        print("t", t, "Calcluated with ", segmentRuns, " segments ", "Upperbound ", upperBound)
        print("Results=", integIter)

    return integIter
