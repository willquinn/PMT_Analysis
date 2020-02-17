import sys
sys.path.insert(1, '..')

import matplotlib.pyplot as plt
import numpy as np
from functions.other_functions import io_parse_arguments, process_date

def main():
    args = io_parse_arguments()
    data_file_Ch0_name = args.i + "resolution_vs_date_Ch0.txt"
    data_file_Ch1_name = args.i + "resolution_vs_date_Ch1.txt"

    try:
        data_file_Ch0 = open(data_file_Ch0_name, 'r')
        data_file_Ch1 = open(data_file_Ch1_name, 'r')
    except FileNotFoundError as fnf_error:
        print(fnf_error)
        raise Exception("Error opening data file")

    data_Ch0 = np.loadtxt(data_file_Ch0_name, delimiter=',', dtype={'names': ('unix', 'exposure', 'date', 'timestamp', 'mu', 'mu_err', 'sigma', 'sigma_err', 'res', 'chi_err'), 'formats': ('i4', 'i4','i4', 'i4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4')}, unpack=True)

    data_Ch1 = np.loadtxt(data_file_Ch1_name, delimiter=',',
                          dtype={'names': ('unix', 'date', 'timestamp', 'mu', 'mu_err', 'sigma', 'sigma_err', 'res', 'chi_err'),
                                 'formats': ('i4', 'i4', 'i4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4')}, unpack=True)

    x_ch0 = process_date(data_Ch0[2])
    x_ch1 = process_date(data_Ch1[1])

    plt.errorbar(x_ch0, data_Ch0[4], data_Ch0[5], fmt='.', color='b')
    plt.title("1MeV PMT Charge Ch0")
    plt.xlabel("Exposure Time /days")
    plt.ylabel("1Mev Charge /pC")
    plt.grid(True)
    plt.show()

    plt.errorbar(x_ch1, data_Ch1[3], data_Ch1[4], fmt='.', color='b')
    plt.title("1MeV PMT Charge Ch1")
    plt.xlabel("Exposure Time /days")
    plt.ylabel("1Mev Charge /pC")
    plt.grid(True)
    plt.show()

    fig, ax1 = plt.subplots()
    ax1.set_title("1MeV PMT Resolution Ch0")
    ax1.grid(True)
    color = 'tab:red'
    plt.errorbar(x_ch0, data_Ch0[8],
                 data_Ch0[8] * np.sqrt((data_Ch0[7] / data_Ch0[6]) ** 2 + (data_Ch0[5] / data_Ch0[4]) ** 2), fmt='.',
                 color=color)
    ax1.set_xlabel("Exposure Time /days")
    ax1.set_ylabel("1Mev Resolution /%", color=color)
    ax1.set_ylim(2, 4.25)
    ax1.tick_params(axis='y', labelcolor=color)

    color = 'tab:blue'
    ax2 = ax1.twinx()
    ax2.plot(x_ch0, data_Ch0[9], '.', color=color)
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

    plt.show()


if __name__ == '__main__':
    main()
