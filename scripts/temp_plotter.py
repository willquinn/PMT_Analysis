import sys
sys.path.insert(1, '..')

import matplotlib.pyplot as plt
import numpy as np
from functions.other_functions import io_parse_arguments, process_date, process_exposure, linear, chi2
from scipy.optimize import curve_fit

def main():
    args = io_parse_arguments()
    data_file_Ch0_name = args.i + "apulseNUM_vs_date_Ch0.txt"
    data_file_Ch1_name = args.i + "apulseNUM_vs_date_Ch1.txt"

    try:
        data_file_Ch0 = open(data_file_Ch0_name, 'r')
        data_file_Ch1 = open(data_file_Ch1_name, 'r')
    except FileNotFoundError as fnf_error:
        print(fnf_error)
        raise Exception("Error opening data file")

    data_Ch0 = np.loadtxt(data_file_Ch0_name, delimiter=',', dtype={'names': ('voltage', 'unix', 'date', 'timestamp', 'perc_l', 'std_l', 'perc_m', 'std_m', 'perc_h', 'std_h'),
                                 'formats': ('i4', 'i4', 'i4', 'i4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4')}, unpack=True)

    data_Ch1 = np.loadtxt(data_file_Ch1_name, delimiter=',', dtype={'names': ('voltage', 'unix', 'date', 'timestamp', 'perc_l', 'std_l', 'perc_m', 'std_m', 'perc_h', 'std_h'),
                                 'formats': ('i4', 'i4', 'i4', 'i4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4')}, unpack=True)

    x_ch0 = process_exposure(data_Ch0[2])
    x_ch1 = process_date(data_Ch1[2])
    x_array = np.linspace(-30, 145, 10000)

    #plt.errorbar(x_ch0, data_Ch0[3], data_Ch0[4], fmt='.', color='b')
    plt.plot(x_ch0, data_Ch0[4], '.', color='r', label='low cut')
    plt.plot(x_ch0, data_Ch0[6], '.', color='g', label='medium cut')
    plt.plot(x_ch0, data_Ch0[8], '.', color='b', label='high cut')
    plt.title("Average number of afterpulses per waveform ch0")
    plt.xlabel("Exposure Time atm-days")
    plt.ylabel("Percentage number of afterpulses /%")
    plt.legend()
    #plt.ylim(20, 60)
    #plt.xscale('log')
    plt.grid(True)
    plt.show()

    # plt.errorbar(x_ch0, data_Ch0[3], data_Ch0[4], fmt='.', color='b')
    plt.plot(x_ch1, data_Ch1[4], '.', color='r', label='low cut')
    plt.plot(x_ch1, data_Ch1[6], '.', color='g', label='medium cut')
    plt.plot(x_ch1, data_Ch1[8], '.', color='b', label='high cut')
    plt.title("Average number of afterpulses per waveform ch1")
    plt.xlabel("Exposure Time /days")
    plt.ylabel("Percentage number of afterpulses /%")
    plt.legend()
    #plt.ylim(20, 60)
    # plt.xscale('log')
    plt.grid(True)
    plt.show()

    '''p_guess = [0.02420914, 3]

    popt_0, pcov_0 = curve_fit(linear, x_ch0, data_Ch0[7],
                               sigma=(data_Ch0[7] * np.sqrt((data_Ch0[6] / data_Ch0[5]) ** 2 + (data_Ch0[4] / data_Ch0[3]) ** 2)),
                               p0=p_guess)
    popt_1, pcov_1 = curve_fit(linear, x_ch1, data_Ch1[7],
                               sigma=(data_Ch1[7] * np.sqrt((data_Ch1[6] / data_Ch1[5]) ** 2 + (data_Ch1[4] / data_Ch1[3]) ** 2)),
                               p0=p_guess)

    chi2_0 = chi2(linear(x_ch0,*popt_0),
                  (data_Ch0[7] * np.sqrt((data_Ch0[6] / data_Ch0[5]) ** 2 + (data_Ch0[4] / data_Ch0[3]) ** 2)),
                   data_Ch0[7],
                  2)
    chi2_1 = chi2(linear(x_ch1, *popt_1),
                  (data_Ch1[7] * np.sqrt((data_Ch1[6] / data_Ch1[5]) ** 2 + (data_Ch1[4] / data_Ch1[3]) ** 2)),
                   data_Ch1[7], 2)

    print(popt_0)
    print(popt_1)

    fig, ax1 = plt.subplots()
    ax1.set_title("1MeV PMT Resolution Ch0")
    ax1.grid(True)
    color = 'tab:red'
    plt.errorbar(x_ch0, data_Ch0[7],
                 data_Ch0[7] * np.sqrt((data_Ch0[6] / data_Ch0[5]) ** 2 + (data_Ch0[4] / data_Ch0[3]) ** 2), fmt='.',
                 color=color)
    #plt.plot(x_array, linear(x_array, *popt_0),'g')
    #ax1.text(-20, 2.1,
    #         'linear fit: y = ({:.4f}±{:.4f})*x + ({:.3f}±{:.3f}) chi2R: {:.3f}'.format(popt_0[0], np.sqrt(pcov_0[0][0]),
    #                                                                                popt_0[1], np.sqrt(pcov_0[1][1]),
    #                                                                                chi2_0))

    ax1.set_xlabel("Exposure Time atm-days")
    ax1.set_ylabel("1Mev Resolution /%", color=color)
    ax1.set_ylim(2, 4.25)
    ax1.tick_params(axis='y', labelcolor=color)

    color = 'tab:blue'
    ax2 = ax1.twinx()
    ax2.plot(x_ch0, data_Ch0[7], '.', color=color)
    ax2.set_ylabel('Chi2', color=color)  # we already handled the x-label with ax1
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.set_ylim(0, 20)
    fig.tight_layout()

    plt.show()

    fig, ax1 = plt.subplots()
    ax1.set_title("1MeV PMT Resolution Ch1")
    ax1.grid(True)
    color = 'tab:red'
    plt.errorbar(x_ch1, data_Ch1[7],
                 data_Ch1[7] * np.sqrt((data_Ch1[6] / data_Ch1[5]) ** 2 + (data_Ch1[4] / data_Ch1[3]) ** 2), fmt='.',
                 color=color)
    plt.plot(x_array, linear(x_array, *popt_1), 'g')
    ax1.text(-20,2.25,
             'linear fit: y = ({:.4f}±{:.4f})*x + ({:.3f}±{:.3f}) chi2R: {:.3f}'.format(popt_1[0], np.sqrt(pcov_1[0][0]),
                                                                                    popt_1[1], np.sqrt(pcov_1[1][1]),
                                                                                    chi2_1))

    ax1.set_xlabel("Exposure Time /days")
    ax1.set_ylabel("1Mev Resolution /%", color=color)
    ax1.set_ylim(2, 4.25)
    ax1.tick_params(axis='y', labelcolor=color)

    color = 'tab:blue'
    ax2 = ax1.twinx()
    ax2.plot(x_ch1, data_Ch1[8], '.', color=color)
    ax2.set_ylabel('Chi2', color=color)  # we already handled the x-label with ax1
    ax2.tick_params(axis='y', labelcolor=color)
    #ax2.set_ylim(0, 20)
    fig.tight_layout()

    plt.show()'''


if __name__ == '__main__':
    main()
