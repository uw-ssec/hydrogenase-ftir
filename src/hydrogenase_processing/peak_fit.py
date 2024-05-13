from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt

def gaussian(x, *params):

    y = np.zeros_like(x)

    for i in range(0, len(params), 3):
        amplitude = params[i]
        center = params[i+1]
        sigma = params[i+2]
        y += amplitude*(1/(sigma* np.sqrt(2*np.pi)))*np.exp((-1/2)*((x - center)/ sigma)**2) 
    return y

def lorentzian(x, *params):
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        amplitude = params[i]
        center = params[i+1]
        sigma = params[i+2]
        y += amplitude*sigma**2/((x-center)**2+sigma**2)
    return y


def peak_fit(fit_function, guess, x_wavenumber, y_absorbance, showplot = True):    
#2091.614748864712, 2080.0497375656423,
    #guess = [0.3, 1939, 2, 0.14, 1958, 2, 0.014, 2080.05, 2, 0.24, 2091.6, 2]
    params, _ = curve_fit(fit_function, x_wavenumber, y_absorbance, guess)

    if showplot:
        plt.plot(x_wavenumber,y_absorbance, label = 'baseline corrected data')
        plt.plot(x_wavenumber, fit_function(x_wavenumber, *params), label = f'{fit_function} fit curve')
        plt.legend()
    return params