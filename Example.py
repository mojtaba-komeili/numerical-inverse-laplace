import numpy as np
import matplotlib.pyplot as plt
from NumInvLaplace import *

# Here we check the inverse transfomr for a known function
# We know that L{ 1 / (1+s) } = e^(-t) is the function that we test
# So, our unit test is finding the inverse of f(x) = 1/(1+x)

tList = np.arange(0.1, 2, 0.1) # predefined range

######### Plotting the analytical inverse
# Simply plotting exp(-t)
yList_analytical = np.exp(-tList)
plt.subplot(2, 1, 1)
plt.plot(tList, yList_analytical)

######### Plotting the numerical inverse
# Plotting the data points that are calculated analytically

f = lambda x: 1 / (1+x) # The 

yList_numerical = []
for t in tList:
    y = inverseLaplace(f, t, 1)
    yList_numerical.append(y)

plt.subplot(2, 1, 2)
plt.plot(tList, yList_numerical)