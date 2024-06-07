from scipy.optimize import curve_fit
from scipy.signal import peak_widths
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


def peak_fit(fit_function, x_wavenumber, y_absorbance, peak_index, showplot = True): 

    get_peak_wid = peak_widths(y_absorbance, peak_index, rel_height=1) 
    all_peak_wid = get_peak_wid[0]
    guess = []
    for i, idx in enumerate(peak_index):
        peak_height = y_absorbance[idx]
        guess.append(peak_height)
        peak_center = x_wavenumber[idx]
        guess.append(peak_center)
        peak_wid = all_peak_wid[i]
        guess.append(peak_wid)

    if fit_function == 'Gaussian':
        objective_function = gaussian
    elif fit_function == 'Lorentzian':
        objective_function = lorentzian
    else:
        raise ValueError("fit_function must be 'Gaussian' or 'Lorentzian'")  # Handle invalid input

    params, _ = curve_fit(objective_function, x_wavenumber, y_absorbance, guess)

    # Calculate RMSE
    predicted_absorbance = objective_function(x_wavenumber, *params)
    residuals = y_absorbance - predicted_absorbance
    rmse = np.sqrt(np.mean(residuals**2))

    if showplot:
        plt.plot(x_wavenumber,y_absorbance, label = 'baseline corrected data')
        plt.plot(x_wavenumber, predicted_absorbance, label = f'{fit_function} curve fit')
        plt.legend()
        
    return params, rmse