from scipy.optimize import curve_fit
import numpy as np


def linear_func_to_fit(x, a, b):
    return b * x + a


def power_law_curve_fit(degree_distribution_list):
    # pass in the degree distribution as a list of [freq. k=0, freq. k=1, ... , freq. k=k_max]
    # then this function will shift it into the (x>0, y>0) interior and fit a power law decay

    # default (failure) values
    found_non_zero = False
    start_x = 0
    a, b = [0, 0]
    opt_para_covariance = [0.0, 0.0]
    fit_success = 0

    # first identify the index of the first non-zero element - we will only fit this (assumed downward) slope
    for x_possible in range(len(degree_distribution_list)):
        if degree_distribution_list[x_possible] > 0:
            start_x = x_possible
            found_non_zero = True
            break

    if found_non_zero:
        # slice y-axis vector to just the relevant part
        non_zero_distribution = degree_distribution_list[start_x:]
        x0 = np.linspace(start_x, start_x + len(non_zero_distribution) - 1, len(non_zero_distribution))
        x1 = x0 + 0.00001  # shift right to address log(degree=0) possibly being counted and having non-zero value
        x2 = np.log(x1)
        # shift the y-values up from zero to be all strictly positive
        y0 = np.array(non_zero_distribution)
        min_y = np.min(y0)
        y1 = y0 + np.abs(min_y) + 0.00001
        y2 = np.log(y1)
        # attempt the linear curve fit of the transformed function
        try:
            optimal_parameters, opt_para_covariance = curve_fit(linear_func_to_fit, x2, y2)
            a, b = [np.exp(optimal_parameters[0]), optimal_parameters[1]]
            fit_success = 1
        except RuntimeError:
            # optimal parameters are not found
            pass

    return fit_success, a, b, opt_para_covariance
