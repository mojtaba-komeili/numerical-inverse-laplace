# Numerical inverse Laplace transformation

# Integration
######################################
import scipy
from scipy import *
from scipy import integrate

# relError_tolerance: The tolerance for accuracy
relError_tolerance = 1e-8

# echo_out_process: A flag
# set it to non-zero values to echo out the mid results.
echo_out_process = 0

def inverseLaplace(func, t, a):
    boundSteps = 5.0 / sqrt(t)
    lowerBound = 0
    upperBound = lowerBound + boundSteps

    integF    = lambda u: real(func(complex(a+1j*u))) * cos(u*t)

    lastInegIter = 0
    integIter    = (2.0*exp(a*t) / pi) * integrate.quad(integF, lowerBound, upperBound)[0]

    segmentRuns = 1
    relError    = abs(abs(integIter) - abs(lastInegIter)) / abs(lastInegIter+1e-15)
    while relError > relError_tolerance:
        segmentRuns  = segmentRuns + 1
        lowerBound   = upperBound
        upperBound   = upperBound + boundSteps
        lastInegIter = integIter
        integIter    = integIter + (2.0*exp(a*t) / pi) * integrate.quad(integF, lowerBound, upperBound)[0]
        relError     = abs(abs(integIter) - abs(lastInegIter)) / (abs(lastInegIter)+1e-15)

    if (echo_out_process == 0):
        print "t", t, "Calcluated with ", segmentRuns, " segments ", "Upperbound ", upperBound
        print "Results=", integIter

    return integIter  
