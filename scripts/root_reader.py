import sys
sys.path.insert(1, '..')

import matplotlib.pyplot as plt
import numpy as np
import ROOT
from scipy.optimize import curve_fit
from functions.other_functions import parse_arguments, get_date_time, get_voltage
from PMT_Array import PMT_Array
from PMT_Waveform import PMT_Waveform


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
    voltage = int(get_voltage(input_data_file_name))
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
        if voltage == 1000:
            pmt_array.set_pmt_templates(
                ["/unix/nemo4/PMT_He_Study_nemo4/Templates/new/191008_A1000_B1000_templates.root"],
                ["A1000_B1000_Ch0_Template", "A1000_B1000_Ch1_Template"])
        elif voltage == 1400:
            pmt_array.set_pmt_templates(
                ["/unix/nemo4/PMT_He_Study_nemo4/Templates/new/190621_A1400_B1400_templates.root"],
                ["A1400_B1400_Ch0_Template", "A1400_B1400_Ch1_Template"])

    i = 0
    for event in tree:
        OM_ID = event.OM_ID
        event_num = event.event_num
        pulse_time = event.pulse_time
        pulse_amplitude = event.pulse_amplitude
        pulse_charge = event.pulse_charge
        pulse_baseline = event.pulse_baseline
        waveform = event.waveform

        pmt_waveform = PMT_Waveform(waveform, pmt_array.get_pmt_object_number(OM_ID))

        if pmt_waveform.get_pulse_trigger():
            pmt_waveform.fill_pmt_hists()

        del pmt_waveform

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
