import sys
sys.path.insert(1, '..')

import matplotlib.pyplot as plt
import numpy as np
import ROOT
from scipy.optimize import curve_fit
from functions.other_functions import parse_arguments, get_date_time
from PMT_Array import PMT_Array
from PMT_Waveform import PMT_Waveform


def chi2(y_obs, y_err, y_exp, n_par):
    chi2 = 0
    ndof = len(y_obs) - n_par - 1
    for i in range(len(y_exp)):
        chi2 += ((y_exp[i] - y_obs[i])/y_err[i])**2
    chi2 = chi2/ndof
    return chi2


def gaus(x, mu, sig):
    return np.exp(-( ( ( x - mu )/sig )**2 )/2)/( np.sqrt( 2*np.pi )*sig )


def fit_0(x, mu, sig, a):
    y = []
    for i in range(len(x)):
        calc =  a * (
                ( 7.08*gaus(x[i], mu, sig ) +
                ( 1.84*gaus(x[i], mu*(1 + 72.144/975.651), sig*1.036) ) +
                ( 0.44*gaus(x[i], mu*(1 + 84.154/975.651), sig*1.042) ) ) +
                0.464 *
                ( np.exp( 0.254*x[i] )/( 1 + np.exp( ( x[i] - 28.43 )/2.14) ) )
        )
        y.append(calc)
    return y


def fit_1(x, mu, sig, a):
    y = []
    for i in range(len(x)):
        calc = a * (
                (7.08 * gaus(x[i], mu, sig) +
                 (1.84 * gaus(x[i], mu * (1 + 72.144 / 975.651), sig * 1.036)) +
                 (0.44 * gaus(x[i], mu * (1 + 84.154 / 975.651), sig * 1.042))) +
                0.515 *
                (np.exp(0.2199 * x[i]) / (1 + np.exp((x[i] - 31.68) / 2.48)))
        )
        y.append(calc)
    return y


def fit(x, mu, sig, a, b, c, d, e):
    y = []
    for i in range(len(x)):
        calc = a * (
                (7.08 * gaus(x[i], mu, sig) +
                 (1.84 * gaus(x[i], mu * (1 + 72.144 / 975.651), sig * 1.036)) +
                 (0.44 * gaus(x[i], mu * (1 + 84.154 / 975.651), sig * 1.042))) +
                b *
                (np.exp(c * x[i]) / (1 + np.exp((x[i] - d) / e)))
        )
        y.append(calc)
    return y


def main():

    # Handle the input arguments:
    ##############################
    args = parse_arguments()
    input_data_file_name = args.i
    config_file_name = args.c
    sweep_bool = args.sweep
    recreate_bool = args.r
    bismuth_bool = args.f
    ##############################

    # Do some string manipulation to get the date and time from the file name
    #################################################
    date, time = get_date_time(input_data_file_name)
    #################################################

    # Check to see if the data file exists
    try:
        print(">>> Reading data from file: {}".format(input_data_file_name))
        date_file = open(input_data_file_name, 'r')
    except FileNotFoundError as fnf_error:
        print(fnf_error)
        raise Exception("Error opening data file {}".format(input_data_file_name))

    file = ROOT.TFile(input_data_file_name, "READ")
    file.cd()

    tree = file.T
    tree.Print()

    topology = [2, 1]
    pmt_array = PMT_Array(topology, date + "_" + time)
    pmt_array.set_pmt_id("GAO607_" + date + "_" + time, 0)
    pmt_array.set_pmt_id("GAO612_" + date + "_" + time, 1)

    # Set the cuts you wish to apply
    # If you don't do this the defaults are used
    if config_file_name is not None:
        pmt_array.apply_setting(config_file_name)

    if sweep_bool == "True":
        pmt_array.set_sweep_bool(True)
        pmt_array.set_pmt_templates("/Users/willquinn/Documents/PhD/PMT_Permeation_Project/data/21.06.19_A1400_B1400_t1130_templates.root", "Template_Waveform_Channel0_A1400_B1400_t1130")

    i = 0
    for event in tree:
        pmt_waveform = PMT_Waveform(event.waveform, pmt_object)

        print(event.OM_ID, event.event_num, event.pulse_time, event.pulse_charge, event.pulse_amplitude, event.pulse_baseline)
        #charge.append(event.pulse_charge)
        #x.append(i)
        #i += 1

    '''y_array_og, bins = np.histogram(charge, bins=150, range=(0,50))
    y_err_og = []
    x_array_og = bins[1:] - bins[1]/2

    range_lower = np.argmax(y_array_og)-int(5/bins[1])
    range_higher = np.argmax(y_array_og)+int(10/bins[1])

    y_array = y_array_og[range_lower:range_higher]
    x_array = x_array_og[range_lower:range_higher]
    for i in range(y_array_og.size):
        if y_array_og[i] == 0:
            y_err_og.append(0)
        else:
            y_err_og.append(np.sqrt(y_array_og[i]))
    y_err_og = np.array(y_err_og)

    y_err = y_err_og[range_lower:range_higher]

    est_mu = np.argmax(y_array_og) * bins[1]
    p_guess = [est_mu, 1, 50, 0.5, 0.2, 28, 2]

    print(p_guess)
    p_bounds = [[est_mu-1, 0, 0, 0, 0, 0, 0], [est_mu+1, 10, 1000, 1000, 1000, 1000, 1000]]

    popt, pcov = curve_fit(fit, x_array, y_array, p0=p_guess, bounds=p_bounds, maxfev=5000)

    print("The optimised fitted parameters are: ", popt)
    print("The covariance matrix is: ", pcov)
    for i in range(len(popt)):
        print("Error on parameter {}: {} is {}".format(i, popt[i], np.sqrt(pcov[i][i])))

    chi_2 = chi2(y_array, y_err, fit(x_array, *popt), 3)
    print("The reduced chi2 is: ", chi_2)

    x = np.linspace(np.min(x_array), np.max(x_array), 10000)

    #print(charge)
    plt.hist(charge, bins=150, range=(0,50), facecolor='b', alpha=0.25)
    plt.errorbar(x_array_og,y_array_og,yerr=y_err_og,fmt='.')
    plt.plot(x, fit(x, *popt))
    plt.text(38,1000,"$\chi^2$ = {}".format(round(chi_2,2)))
    plt.grid(True)
    plt.ylabel("Counts")
    plt.xlabel("Charge /pC")
    plt.yscale('log')
    plt.show()'''


if __name__ == '__main__':
    main()
