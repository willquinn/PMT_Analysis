import sys

sys.path.insert(1, '..')

# import ROOT and bash commands
import ROOT
import os
import tqdm

# import python plotting and numpy modules
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

# import stats module
from scipy.optimize import curve_fit

# import custom made classes
from functions.other_functions import sncd_parse_arguments, fit, chi2, process_date, linear
from scr.PMT_Array import PMT_Array


def print_settings(pmt_array: PMT_Array):
    pass


def get_error_A_divide_B(dA, A, dB, B):
    C = A/B
    dC = C * np.sqrt( (dA/A)**2 + (dB/B)**2 )
    return dC


def get_resolution(mu: float, mu_err: float, sig: float, sig_err: float):
    res = 0
    res_err = 0
    if mu == 0:
        pass
    else:
        res = sig/mu
        res_err = get_error_A_divide_B(sig_err, sig, mu_err, mu)
    return res, res_err


def read_tree(date: str, root_file_name: str, pmt_array: PMT_Array, output_file_location: str, output_file_name: str):

    file = ROOT.TFile(root_file_name, "READ")
    file.cd()

    # date = root_file_name.split("_")[0]

    tree = file.T

    # Create the histograms that we will want to store
    # these will be used to extract the resolution
    nbins = pmt_array.get_pmt_object_number(0).get_setting("nbins")

    charge_hist_ch0 = ROOT.TH1D(date + "_" + pmt_array.get_pmt_object_number(0).get_pmt_id() + "_charge_spectrum",
                                date + "_" + pmt_array.get_pmt_object_number(0).get_pmt_id() + "_charge_spectrum",
                                nbins, 0, 60)
    charge_hist_ch1 = ROOT.TH1D(date + "_" + pmt_array.get_pmt_object_number(1).get_pmt_id() + "_charge_spectrum",
                                date + "_" + pmt_array.get_pmt_object_number(1).get_pmt_id() + "_charge_spectrum",
                                nbins, 0, 60)

    charge_hists = ROOT.TList()
    charge_hists.Add(charge_hist_ch0)
    charge_hists.Add(charge_hist_ch1)

    # The amplitudes and baselines are less important
    # but it will be interesting to see have they vary in time
    raw_amplitudes = [[] for i in range(pmt_array.get_pmt_total_number())]
    raw_baselines = [[] for i in range(pmt_array.get_pmt_total_number())]
    raw_charges = [[] for i in range(pmt_array.get_pmt_total_number())]

    fit_parameters = []

    for event in tree:

        # Access the information inside the NTuple
        OM_ID = event.OM_ID
        # event_num = event.event_num
        # pulse_time = event.pulse_time
        pulse_amplitude = event.pulse_amplitude
        pulse_charge = event.pulse_charge
        pulse_baseline = event.pulse_baseline
        apulse_num = event.apulse_num
        apulse_times = np.array(event.apulse_times)
        apulse_amplitudes = np.array(event.apulse_amplitudes)
        apulse_shapes = np.array(event.apulse_shapes)

        # Apply the cuts
        if pmt_array.get_pmt_object_number(OM_ID).get_setting("charge_cut") < pulse_charge:
            continue

        if pulse_amplitude > 0 and pulse_baseline > 0:
            raw_amplitudes[OM_ID].append(pulse_amplitude)
            raw_baselines[OM_ID].append(pulse_baseline)
            raw_charges[OM_ID].append(pulse_charge)

            charge_hists[OM_ID].Fill(pulse_charge)

            filer_list = []
            for i_apulse in range(apulse_num):
                if apulse_shapes[i_apulse] > pmt_array.get_pmt_object_number(OM_ID).get_setting("mf_shape_threshold")\
                        and apulse_amplitudes[i_apulse] > pmt_array.get_pmt_object_number(OM_ID).get_setting("mf_amp_threshold")\
                        and apulse_times[i_apulse] > pmt_array.get_pmt_object_number(OM_ID).get_setting("sweep_range")[0]:
                    filer_list.append(True)
                else:
                    filer_list.append(False)

    # Save a png of the amplitude and baseline for each day
    # Now fit a Bi-207 function to the area spectrum
    my_amp_bins = [i*1000/nbins for i in range(nbins)]
    my_charge_bins = [i*60 / nbins for i in range(nbins)]
    my_baseline_bins = [900 + (i*100/nbins) for i in range(nbins)]
    for i_om in range(pmt_array.get_pmt_total_number()):

        if len(raw_amplitudes[i_om]) < 1:
            continue

        plt.hist(raw_amplitudes[i_om], my_amp_bins, color='g')
        plt.xlabel("Amplitude /mV")
        plt.ylabel("counts")
        plt.title(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() + "_amplitude_spectrum")
        plt.grid()
        plt.yscale('log')
        plt.savefig(output_file_location + "/plots/" + date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() + "_amplitude_spectrum")
        plt.close()

        plt.hist(raw_baselines[i_om], my_baseline_bins, color='r')
        plt.xlabel("Baseline /mV")
        plt.ylabel("counts")
        plt.title(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() + "_baseline")
        plt.grid()
        plt.yscale('log')
        plt.savefig(
            output_file_location + "/plots/" + date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() + "_baseline")
        plt.close()

        bin_contents, bin_edges = np.histogram(raw_charges[i_om], my_charge_bins)
        peak = np.argmax(bin_contents)
        bin_centers = bin_edges[:-1] + np.diff(bin_edges) / 2

        upper_limit = peak + 5 * 60 / nbins
        lower_limit = peak - 3 * 60 / nbins

        p_guess = [bin_centers[peak], 0.04*bin_centers[peak], bin_contents[peak], 0.2, bin_centers[peak] - 5, 2]
        p_bounds = [[bin_centers[peak] - 2, 0.04*bin_centers[peak] - 2*0.04*bin_centers[peak], bin_contents[peak] - 50, 0.3, bin_centers[peak] - 5 - 2, 1],
                    [bin_centers[peak] + 2, 0.04*bin_centers[peak] + 2*0.04*bin_centers[peak], bin_contents[peak] + 50, 0.1, bin_centers[peak] - 5 + 2, 3]
                    ]

        popt, pcov = curve_fit(fit, bin_centers[lower_limit:upper_limit],
                               bin_contents[lower_limit:upper_limit],
                               sigma=np.sqrt(bin_contents[lower_limit:upper_limit]),
                               p0=p_guess, bounds=p_bounds, maxfev=1000)

        chi_2 = chi2(fit(bin_centers[lower_limit:upper_limit], *popt),
                     np.sqrt(bin_contents[lower_limit:upper_limit]),
                     bin_contents[lower_limit:upper_limit],
                     len(popt))

        pars = {
            "mu": popt[0],
            "sig": popt[1],
            "mu_err": np.sqrt(pcov[0,0]),
            "sig_err": np.sqrt(pcov[1,1]),
            "chi2": chi_2
        }
        fit_parameters.append(pars)

        plt.hist(raw_charges, my_charge_bins)
        plt.plot(bin_centers[lower_limit:upper_limit],
                 fit(bin_centers[lower_limit:upper_limit], *popt))
        plt.xlabel("charge /pC")
        plt.ylabel("counts")
        plt.title(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() + "_charge_spectrum")
        plt.grid()
        plt.yscale('log')
        plt.savefig(output_file_location + "/plots/" + date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() + "_charge_spectrum")
        plt.close()

    file.Close()

    output_file = ROOT.TFile(output_file_location + "/ROOT_files/" + output_file_name, "RECREATE")
    output_file.cd()

    charge_hist_ch0.Write()
    charge_hist_ch1.Write()

    output_file.Close()

    return fit_parameters


def main():
    # Handle the input arguments:
    ##############################
    args = sncd_parse_arguments()
    filenames_txt = args.i
    config_file_name = args.c
    output_directory = args.o
    ##############################

    # Do some string manipulation to get the date and time from the file name
    #################################################
    voltage = 1000
    #################################################

    try:
        print(">>> Reading data from file: {}".format(filenames_txt))
        date_file = open(filenames_txt, 'r')
    except FileNotFoundError as fnf_error:
        print(fnf_error)
        raise Exception("Error opening data file {}".format(filenames_txt))

    filenames = np.loadtxt(filenames_txt, delimiter=',', dtype={
        'names': ['filename'],
        'formats': ['S1']}, unpack=True)

    topology = [2, 1]
    pmt_array = PMT_Array(topology, "summary")
    pmt_array.set_pmt_id("GAO607", 0)
    pmt_array.set_pmt_id("GAO612", 1)

    # Set the cuts you wish to apply
    # If you don't do this the defaults are used
    if config_file_name is not None:
        pmt_array.apply_setting(config_file_name)
        print_settings(pmt_array)

    # Set up the containers for the summary
    resolutions = [[] for i in range(pmt_array.get_pmt_total_number())]
    resolutions_err = [[] for i in range(pmt_array.get_pmt_total_number())]
    dates = [[] for i in range(pmt_array.get_pmt_total_number())]
    gains = [[] for i in range(pmt_array.get_pmt_total_number())]
    gains_err = [[] for i in range(pmt_array.get_pmt_total_number())]
    fit_chi2 = [[] for i in range(pmt_array.get_pmt_total_number())]

    for i_file, file in tqdm.tqdm(enumerate(filenames.size)):
        date = file.split("_")[0]

        fit_parameters = read_tree(date, file, pmt_array, "", file.split(".root")[0] + "_output.root")

        for i_om in range(pmt_array.get_pmt_total_number()):

            mu = fit_parameters[i_om]["mu"]
            mu_err = fit_parameters[i_om]["mu_err"]
            sig = fit_parameters[i_om]["sig"]
            sig_err = fit_parameters[i_om]["sig_err"]
            chi_2 = fit_parameters[i_om]["chi2"]

            res, res_err = get_resolution(mu, mu_err, sig, sig_err)

            resolutions[i_om].append(res)
            resolutions_err[i_om].append(res_err)
            dates[i_om].append(date)
            gains.append(mu)
            gains_err[i_om].append(mu_err)
            fit_chi2[i_om].append(chi_2)

        break

    # Plot individual summaries
    for i_om in range(pmt_array.get_pmt_total_number()):

        date = np.array(process_date(dates[i_om]))

        start = np.where(date == 0)
        mid = np.where(date == 98)

        popt, pcov = curve_fit(linear, date[start:], resolutions[i_om][start:],
                               sigma=resolutions_err[i_om][start:], p0=[1, 3], bounds=[[0, 0], [10, 10]])
        x_array = np.linspace(date[start], np.amax(date), 2)
        chi_2 = chi2(resolutions[i_om][start:], resolutions_err[i_om][start:], linear(date[start:], *popt), 2)

        plt.subplot(2, 1, 1)
        plt.errorbar(date[:start + 1], resolutions[i_om][:start + 1], yerr=resolutions_err[i_om][:start + 1],
                     fmt="g.", label="Atmospheric He")
        plt.errorbar(date[start+1:mid + 1], resolutions[i_om][start+1:mid + 1], yerr=resolutions_err[i_om][start+1:mid + 1],
                     fmt="b.", label="1% He")
        plt.errorbar(date[mid+1:], resolutions[i_om][mid+1:], yerr=resolutions_err[i_om][mid+1:],
                     fmt="r.", label="10% He")
        plt.plot(x_array, linear(x_array, *popt), 'k-',
                 label="$y = (${:.1e}$ ± ${:.0e})$\\times x + (${:.1e}$ ± ${:.0e}$) \chi^2_R = {:.2}$".format(popt[0], np.sqrt(pcov[0, 0]),popt[1], np.sqrt(pcov[1, 1]), chi_2))
        plt.xlabel("exposure days relative to 190611")
        plt.ylabel("Resolution at 1MeV /% $\sigma/\mu$")
        plt.title(pmt_array.get_pmt_object_number(i_om).get_pmt_id() + " Resolution vs exposure time")
        plt.axvline(date[start], 0, 100, ls='--', color='k')
        plt.axvline(date[mid], 0, 100, ls='--', color='k')
        plt.grid()
        plt.ylim(2, 4.5)
        plt.xlim(np.amin(date), np.amax(date))
        plt.legend(loc='lower right')

        plt.subplot(2, 1, 2)
        plt.plot(date[:start + 1], fit_chi2[i_om][:start + 1], "g.", label="Atmospheric He")
        plt.plot(date[start + 1:mid + 1], fit_chi2[i_om][start + 1:mid + 1], "b.", label="1% He")
        plt.plot(date[mid + 1:], fit_chi2[i_om][mid + 1:], "r.", label="10% He")
        plt.xlabel("exposure days relative to 190611")
        plt.ylabel("$\chi^2_R$")
        plt.title(pmt_array.get_pmt_object_number(i_om).get_pmt_id() + " Resolution fit $\chi^2_R$ vs exposure time")
        plt.grid()
        plt.axvline(date[start], 0, 100, ls='--', color='k')
        plt.axvline(date[mid], 0, 100, ls='--', color='k')
        plt.legend(loc='lower_right')
        plt.xlim(np.amin(date), np.amax(date))
        plt.ylim(0, 10)

        plt.savefig(output_directory + "/Summary_plots/" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() +
                    "_resolution_vs_time")
        plt.close()

    # Plot ratio
    x_date = []
    ratio = []
    ratio_err = []
    for i in range(len(dates[0])):
        for j in range(len(dates[1])):
            if dates[0][i] == dates[1][j]:
                x_date.append(dates[0][i])
                ratio.append(resolutions[0][i] / resolutions[1][j])
                ratio_err.append(resolutions[0][i] / resolutions[1][j] * np.sqrt(
                    (resolutions_err[0][i] / resolutions[0][i]) ** 2 + (resolutions_err[1][j] / resolutions[1][j]) ** 2))
                break

    popt, pcov = curve_fit(linear, x_date, ratio,
                           sigma=ratio_err, p0=[1, 1], bounds=[[0, 0], [10, 10]])
    x_array = np.linspace(np.amin(x_date), np.amax(x_date), 2)
    chi_2 = chi2(ratio, ratio_err, linear(x_date, *popt), 2)

    plt.errorbar(x_date, ratio, yerr=ratio_err, fmt="k.")
    plt.plot(x_array, linear(x_array, *popt), "g-",
             label="$y = (${:.1e}$ ± ${:.0e})$\\times x + (${:.1e}$ ± ${:.0e}$) \chi^2_R = {:.2}$".format(popt[0], np.sqrt(pcov[0, 0]),popt[1], np.sqrt(pcov[1, 1]), chi_2))
    plt.axvline(98, color="r", ls="--")
    plt.axvline(0, color="b", ls="--")
    plt.xlabel("exposure days relative to 190611")
    plt.ylabel("Ratio res_Ch0/res_Ch1")
    plt.title("Ratio of resolution of CH 0 & 1 vs time")
    plt.grid()
    plt.xlim(np.amin(np.array(x_date)), np.amax(np.array(x_date)))
    plt.ylim(0, 2)
    plt.savefig(output_directory + "/Summary_plots/resolution_ratio_vs_time")
    plt.close()


if __name__ == '__main__':
    main()
