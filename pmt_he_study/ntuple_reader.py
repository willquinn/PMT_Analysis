import sys

sys.path.insert(1, '..')

# import ROOT and bash commands
import ROOT
import tqdm

# import python plotting and numpy modules
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

# import stats module
from scipy.optimize import curve_fit

# import custom made classes
from functions.other_functions import sncd_parse_arguments, fit, chi2, process_date, linear, gaus
from scr.PMT_Array import PMT_Array


def print_settings(pmt_array: PMT_Array):
    pass


def temp_fit(x, a, b, c, d, e):
    y = []
    for i in range(len(x)):
        calc = a* b * (np.exp(c * x[i]) / (1 + np.exp((x[i] - d) / e)))
        y.append(calc)
    return y


def temp_fit1(x, mu, sig, a):
    y = []
    for i in range(len(x)):
        calc = a * (
                (7.08 * gaus(x[i], mu, sig) +
                 (1.84 * gaus(x[i], mu * (1 + 72.144 / 975.651), sig * 1.036)) +
                 (0.44 * gaus(x[i], mu * (1 + 84.154 / 975.651), sig * 1.042))))
        y.append(calc)
    return y


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


def read_tree(root_file_name: str, pmt_array: PMT_Array, output_file_location: str, output_file_name: str):

    file = ROOT.TFile(root_file_name, "READ")
    file.cd()

    date = root_file_name.split("_")[0]
    voltage = int(root_file_name.split("_")[2].split("A")[1])

    if voltage == 1000:
        max_amp = 400
    else:
        max_amp = 1000
    max_charge = 60

    tree = file.T

    # Create the histograms that we will want to store
    # these will be used to extract the resolution
    nbins = pmt_array.get_pmt_object_number(0).get_setting("nbins")

    charge_hist_ch0 = ROOT.TH1D(date + "_" + pmt_array.get_pmt_object_number(0).get_pmt_id() + "_charge_spectrum_" + str(voltage) + "V",
                                date + "_" + pmt_array.get_pmt_object_number(0).get_pmt_id() + "_charge_spectrum_" + str(voltage) + "V",
                                nbins, 0, max_charge)
    charge_hist_ch1 = ROOT.TH1D(date + "_" + pmt_array.get_pmt_object_number(1).get_pmt_id() + "_charge_spectrum_" + str(voltage) + "V",
                                date + "_" + pmt_array.get_pmt_object_number(1).get_pmt_id() + "_charge_spectrum_" + str(voltage) + "V",
                                nbins, 0, max_charge)
    amp_hist_ch0 = ROOT.TH1D(date + "_" + pmt_array.get_pmt_object_number(0).get_pmt_id() + "_amplitude_spectrum_" + str(voltage) + "V",
                             date + "_" + pmt_array.get_pmt_object_number(0).get_pmt_id() + "_amplitude_spectrum_" + str(voltage) + "V",
                             nbins, 0, max_amp)
    amp_hist_ch1 = ROOT.TH1D(date + "_" + pmt_array.get_pmt_object_number(1).get_pmt_id() + "_amplitude_spectrum_" + str(voltage) + "V",
                             date + "_" + pmt_array.get_pmt_object_number(1).get_pmt_id() + "_amplitude_spectrum_" + str(voltage) + "V",
                             nbins, 0, max_amp)
    baseline_ch0 = ROOT.TH1D(date + "_" + pmt_array.get_pmt_object_number(0).get_pmt_id() + "_baseline_distribution_" + str(voltage) + "V",
                             date + "_" + pmt_array.get_pmt_object_number(0).get_pmt_id() + "_baseline_distribution_" + str(voltage) + "V",
                             nbins, 978, 981)
    baseline_ch1 = ROOT.TH1D(date + "_" + pmt_array.get_pmt_object_number(1).get_pmt_id() + "_baseline_distribution_" + str(voltage) + "V",
                             date + "_" + pmt_array.get_pmt_object_number(1).get_pmt_id() + "_baseline_distribution_" + str(voltage) + "V",
                             nbins, 978, 981)

    charge_hists = [charge_hist_ch0, charge_hist_ch1]
    amp_hists = [amp_hist_ch0, amp_hist_ch1]
    baselines = [baseline_ch0, baseline_ch1]

    '''# The amplitudes and baselines are less important
    # but it will be interesting to see have they vary in time
    raw_amplitudes = [[] for i in range(pmt_array.get_pmt_total_number())]
    raw_baselines = [[] for i in range(pmt_array.get_pmt_total_number())]
    raw_charges = [[] for i in range(pmt_array.get_pmt_total_number())]'''

    for event in tree:

        # Access the information inside the NTuple
        OM_ID = event.OM_ID
        if OM_ID == 0 or OM_ID == 1:
            pass
        else:
            continue

        pulse_amplitude = int(event.pulse_amplitude)
        pulse_charge = event.pulse_charge
        pulse_baseline = event.pulse_baseline
        apulse_num = event.apulse_num
        apulse_times = np.array(event.apulse_times)
        apulse_amplitudes = np.array(event.apulse_amplitudes)
        apulse_shapes = np.array(event.apulse_shapes)

        '''raw_amplitudes[OM_ID].append(pulse_amplitude)
        raw_baselines[OM_ID].append(pulse_baseline)
        raw_charges[OM_ID].append(pulse_charge)'''

        charge_hists[OM_ID].Fill(pulse_charge)
        amp_hists[OM_ID].Fill(pulse_amplitude)
        baselines[OM_ID].Fill(pulse_baseline)

        # Now apply new amplitude and shape cuts
        filer_list = []
        for i_apulse in range(apulse_num):
            if apulse_shapes[i_apulse] > pmt_array.get_pmt_object_number(OM_ID).get_setting("mf_shape_threshold")\
                    and apulse_amplitudes[i_apulse] > pmt_array.get_pmt_object_number(OM_ID).get_setting("mf_amp_threshold")\
                    and apulse_times[i_apulse] > pmt_array.get_pmt_object_number(OM_ID).get_setting("sweep_range")[0]:
                filer_list.append(True)
            else:
                filer_list.append(False)

    '''# Save a png of the amplitude and baseline for each day
    # Now fit a Bi-207 function to the area spectrum
    my_amp_bins = [i*max_amp/nbins for i in range(nbins)]
    my_charge_bins = [i*60 / nbins for i in range(nbins)]
    for i_om in range(pmt_array.get_pmt_total_number()):
        if len(raw_amplitudes[i_om]) < 1:
            fit_parameters.append([])
            continue
        plt.hist(raw_amplitudes[i_om], bins=my_amp_bins, color='g')
        plt.xlabel("Amplitude /mV")
        plt.ylabel("counts")
        plt.title(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() + "_amplitude_spectrum")
        plt.grid()
        plt.xlim(0,400)
        #plt.yscale('log')
        plt.savefig(output_file_location + "/plots/" + date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() + "_amplitude_spectrum")
        plt.close()

        plt.hist(raw_baselines[i_om], bins=nbins, color='r')
        plt.xlabel("Baseline /mV")
        plt.ylabel("counts")
        plt.title(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() + "_baseline")
        plt.grid()
        #plt.yscale('log')
        plt.savefig(
            output_file_location + "/plots/" + date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() + "_baseline")
        plt.close()

        bin_contents, bin_edges = np.histogram(raw_charges[i_om], bins=my_charge_bins)
        peak = np.argmax(bin_contents)
        bin_centers = bin_edges[:-1] + np.diff(bin_edges) / 2

        upper_limit = bin_centers[peak] + 5
        lower_limit = bin_centers[peak] - 3

        u_bin = peak + int(10*nbins/60)
        l_bin = peak - int(6*nbins/60)

        p_guess = [bin_centers[peak], 0.04*bin_centers[peak], bin_contents[peak], 0.5, 0.2, bin_centers[peak] - 5, 2]
        # print("p_guess:", p_guess)
        p_bounds = [[p_guess[0] - 2, p_guess[1] - 2 * 0.04 * bin_centers[peak], 0, 0, 0, bin_centers[peak] - 5 - 10, 0],
                    [p_guess[0] + 2, p_guess[1] + 2 * 0.04 * bin_centers[peak], 10000, 1, 1, bin_centers[peak] - 5 + 10, 50]
                    ]
        # print("p_bounds:", p_bounds)

        print("x:", bin_centers[l_bin:u_bin])
        print("y:", bin_contents[l_bin:u_bin])
        print("err:", np.sqrt(bin_contents[l_bin:u_bin])

        try:
            popt, pcov = curve_fit(fit, bin_centers[l_bin:u_bin],
                                   bin_contents[l_bin:u_bin],
                                   sigma=np.sqrt(bin_contents[l_bin:u_bin]),
                                   p0=p_guess, bounds=p_bounds, maxfev=1000000)
        except:
            fit_parameters.append([])
            continue
        # print("popt:", popt)

        chi_2 = chi2(fit(bin_centers[l_bin:u_bin], *popt),
                     np.sqrt(bin_contents[l_bin:u_bin]),
                     bin_contents[l_bin:u_bin],
                     len(popt))

        pars = {
            "mu": popt[0],
            "sig": popt[1],
            "mu_err": np.sqrt(pcov[0,0]),
            "sig_err": np.sqrt(pcov[1,1]),
            "chi2": chi_2
        }
        fit_parameters.append(pars)

        fig, ax = plt.subplots()
        plt.hist(raw_charges[i_om], bins=my_charge_bins, color='g', alpha=0.5)
        plt.plot(bin_centers[l_bin:u_bin],
                 fit(bin_centers[l_bin:u_bin], *popt), 'r-', label='CE: ~1MeV Peak Fit + BKGD')
        plt.plot(bin_centers[l_bin-10:u_bin+10],
                 temp_fit(bin_centers, popt[2], popt[3], popt[4], popt[5], popt[6])[l_bin-10:u_bin+10], 'b,-', label='BKGD')
        plt.plot(bin_centers,
                 temp_fit1(bin_centers, popt[0], popt[1], popt[2]), 'p-')
        plt.xlabel("charge /pC")
        plt.ylabel("counts")
        plt.title(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() + "_charge_spectrum")
        plt.grid()
        textstr = '\n'.join((r'$\mu=%.2f ± %.2f$' % (pars['mu'], pars['mu_err']),
                             r'$\sigma=%.2f ± %.2f$' % (pars['sig'], pars['sig_err']),
                             r'$\mathrm{Amp}=%.2f$' % (popt[2], ),
                             r'$\chi^2_R=%.2f$' % (chi_2,)
                             ))
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        plt.text(0.75, 0.95, textstr, transform=ax.transAxes, fontsize=10, verticalalignment='top', bbox=props)
        plt.yscale('log')
        plt.ylim(1, 10000)
        plt.xlim(0, 60)
        plt.legend(loc='lower right')
        plt.savefig(output_file_location + "/plots/" + date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() + "_charge_spectrum")
        plt.close()'''

    output_file = ROOT.TFile(output_file_location + "/ROOT_files/" + str(voltage) + "V/" + output_file_name, "RECREATE")
    output_file.cd()

    for i_om in range(pmt_array.get_pmt_total_number()):
        if charge_hists[i_om].GetEntries() > 0:
            charge_hists[i_om].Write()
            amp_hists[i_om].Write()
            baselines[i_om].Write()
        else:
            pass

    file.Close()
    output_file.Close()


def main():
    # Handle the input arguments:
    ##############################
    args = sncd_parse_arguments()
    input_directory = args.i
    config_file_name = args.c
    output_directory = args.o
    summary_bool = args.s
    ##############################

    filenames_txt = input_directory + "/1000V_filenames.txt"

    try:
        print(">>> Reading data from file: {}".format(filenames_txt))
        date_file = open(filenames_txt, 'r')
    except FileNotFoundError as fnf_error:
        print(fnf_error)
        raise Exception("Error opening data file {}".format(filenames_txt))

    filenames = np.loadtxt(filenames_txt, delimiter=',', dtype={
        'names': ['filename'],
        'formats': ['S100']}, unpack=True)

    topology = [2, 1]
    pmt_array = PMT_Array(topology, "summary")
    pmt_array.set_pmt_id("GAO607", 0)
    pmt_array.set_pmt_id("GAO612", 1)

    # Set the cuts you wish to apply
    # If you don't do this the defaults are used
    if config_file_name is not None:
        pmt_array.apply_setting(config_file_name)
        print_settings(pmt_array)

    '''# Set up the containers for the summary
    resolutions = [[] for i in range(pmt_array.get_pmt_total_number())]
    resolutions_err = [[] for i in range(pmt_array.get_pmt_total_number())]
    dates = [[] for i in range(pmt_array.get_pmt_total_number())]
    gains = [[] for i in range(pmt_array.get_pmt_total_number())]
    gains_err = [[] for i in range(pmt_array.get_pmt_total_number())]
    fit_chi2 = [[] for i in range(pmt_array.get_pmt_total_number())]'''

    for i_file in tqdm.tqdm(range(filenames.size)):
        file = filenames[i_file][0].decode("utf-8")

        read_tree(input_directory + "/" + file, pmt_array, output_directory, file.split(".root")[0] + "_output.root")

        '''for i_om in range(pmt_array.get_pmt_total_number()):

            if len(fit_parameters[i_om]) > 0:
                pass
            else:
                continue

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
            fit_chi2[i_om].append(chi_2)'''


    '''# Plot individual summaries
    for i_om in range(pmt_array.get_pmt_total_number()):

        if len(fit_parameters[i_om]) > 0:
            pass
        else:
            continue
        print(dates[i_om])
        date = process_date(np.array(dates[i_om]))
        print(date)

        start = np.where(date == 0)[0]
        print(start)
        mid = np.where(date == 98)[0]
        print(mid)

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
    plt.close()'''


if __name__ == '__main__':
    main()
