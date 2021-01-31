#constraining initial core mass distribution and atmospheric mass fraction distribution
import numpy as np
import matplotlib.pyplot as plt
import warnings
from astropy import constants as const
from scipy.optimize import minimize
import warnings
from scipy.optimize import fsolve

#aim
#draw core mass sample from a gaussian
#Uniformly distributed logR logP plane.
#constrain core mass and atmospheric fraction distribution to ensure : (1) sub saturn desert is present (2) current mass-radius relationship is matched 
#(3) xi_sqr minimization with occurence rate


