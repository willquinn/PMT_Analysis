import sys
sys.path.insert(1, '..')

import matplotlib.pyplot as plt
import numpy as np
from functions.other_functions import io_parse_arguments, process_date, process_exposure, linear, chi2
from scipy.optimize import curve_fit


def my_func(t, pe, n, c):
    return pe*c*t/( ( (1/(c*t)) + 1 ) **n )

def my_func_der(t, pe, n, c):
    return pe*c/( ( (1/(c*t)) + 1 ) **n ) + pe*n/(t * ( ( (1/(c*t)) + 1 ) **n ))

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

    apulse_data_Ch0 = np.loadtxt(data_file_Ch0_name, delimiter=',', dtype={'names': ('voltage', 'unix', 'date', 'timestamp', 'perc_l', 'std_l', 'perc_m', 'std_m', 'perc_h', 'std_h'),
                                 'formats': ('i4', 'i4', 'i4', 'i4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4')}, unpack=True)

    apulse_data_Ch1 = np.loadtxt(data_file_Ch1_name, delimiter=',', dtype={'names': ('voltage', 'unix', 'date', 'timestamp', 'perc_l', 'std_l', 'perc_m', 'std_m', 'perc_h', 'std_h'),
                                 'formats': ('i4', 'i4', 'i4', 'i4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4')}, unpack=True)

    print(np.sort(apulse_data_Ch0[2]))
    x_ch0 = process_date(apulse_data_Ch0[2])
    x_ch1 = process_date(apulse_data_Ch1[2])
    x_array = np.linspace(-30, 145, 10000)

    #plt.errorbar(x_ch0, apulse_data_Ch0[3], apulse_data_Ch0[4], fmt='.', color='b')
    plt.plot(x_ch0, apulse_data_Ch0[4], '.', color='r', label='low cut')
    plt.plot(x_ch0, apulse_data_Ch0[6], '.', color='g', label='medium cut')
    plt.plot(x_ch0, apulse_data_Ch0[8], '.', color='b', label='high cut')
    plt.title("Average number of afterpulses per waveform ch0")
    plt.xlabel("Exposure Time days")
    plt.ylabel("Percentage number of afterpulses /%")
    #plt.legend()
    #plt.ylim(20, 60)
    #plt.xscale('log')
    plt.grid(True)
    plt.savefig("/Users/willquinn/Desktop/apulse_vs_date_time_ch0")
    plt.show()

    # plt.errorbar(x_ch0, apulse_data_Ch0[3], apulse_data_Ch0[4], fmt='.', color='b')
    plt.plot(x_ch1, apulse_data_Ch1[4], '.', color='r', label='low cut')
    plt.plot(x_ch1, apulse_data_Ch1[6], '.', color='g', label='medium cut')
    plt.plot(x_ch1, apulse_data_Ch1[8], '.', color='b', label='high cut')
    plt.title("Average number of afterpulses per waveform ch1")
    plt.xlabel("Exposure Time /days")
    plt.ylabel("Percentage number of afterpulses /%")
    #plt.legend()
    #plt.ylim(20, 60)
    # plt.xscale('log')
    plt.grid(True)
    plt.savefig("/Users/willquinn/Desktop/apulse_vs_date_time_ch1")
    plt.show()

    '''t = np.linspace(100,13000000,13000000)
    pi = my_func(t, 10000, 100, 0.0000001)
    dtpi = my_func_der(t, 10000, 100, 0.0000001)
    plt.plot(t,pi,'b')
    plt.yscale('log')
    plt.xlabel("Time /s")
    plt.ylabel("$p_i$ /Pa")
    plt.grid(True)
    plt.show()
    plt.plot(t, pi, 'b')
    plt.xlabel("Time /s")
    plt.ylabel("$p_i$ /Pa")
    plt.grid(True)
    plt.show()
    plt.plot(t,dtpi, 'r')
    plt.yscale('log')
    plt.xlabel("Time /s")
    plt.ylabel("$\\frac{\delta p_i}{\delta t}$ /Pa")
    plt.show()'''


if __name__ == '__main__':
    main()
