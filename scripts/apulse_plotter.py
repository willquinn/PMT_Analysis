import sys
sys.path.insert(1, '..')

import numpy as np
from functions.data_reader_functions import read_filenames
from functions.other_functions import io_parse_arguments, process_date, process_exposure
import matplotlib.pyplot as plt
import ROOT


def main():
    args = io_parse_arguments()
    input_data_filename = args.i
    data_file_list = read_filenames(input_data_filename)
    PATH = input_data_filename.split("filenames")[0]

    print(input_data_filename)

    date_list_ch0 = []
    date_list_ch1 = []
    num_apulse_ch0_list = []
    num_apulse_ch1_list = []

    ch0_num_apulse = []
    ch0_num_apulse_err = []
    ch1_num_apulse = []
    ch1_num_apulse_err = []

    for data_file_index, data_file_name in enumerate(data_file_list):
        try:
            _file = open(PATH + data_file_name, 'r')
        except:
            continue

        file = ROOT.TFile(PATH + data_file_name, "READ")
        d0 = file.GetDirectory("0")
        d1 = file.GetDirectory("1")

        hist_ch0 = d0.Get("h_afterpulse_times_Ch0;1")
        ref_hist0 = d0.Get("h_num_afterpulses_Ch0;1")
        hist_ch1 = d1.Get("h_afterpulse_times_Ch1;1")
        ref_hist1 = d1.Get("h_num_afterpulses_Ch1;1")

        try:
            hist_ch0.GetEntries()
            hist_ch1.GetEntries()
            print(data_file_index)
            print("CH0 Total number of afterpulse: ", hist_ch0.GetEntries())
            print("CH0 Total number of waveforms: ", hist_ch0.GetEntries())
            print("CH1 Total number of afterpulse: ", hist_ch1.GetEntries())
            print("CH1 Total number of waveforms: ", ref_hist1.GetEntries())

        except:
            continue

        temp_ch0 = 0
        temp_ch1 = 0
        date = int(data_file_name.split('_')[0])

        # Count the bins you want to include
        for i_bin in range(hist_ch0.GetNbinsX()):
            if 1500 < i_bin*hist_ch0.GetBinWidth(i_bin) <= 2600:
                temp_ch0 += hist_ch0.GetBinContent(i_bin)
        for i_bin in range(hist_ch1.GetNbinsX()):
            if 1500 < i_bin*hist_ch1.GetBinWidth(i_bin) <= 2600:
                temp_ch1 += hist_ch1.GetBinContent(i_bin)

        if hist_ch0.GetEntries() > 0:
            num_apulse_ch0_list.append(100 * temp_ch0 / hist_ch0.GetEntries())
            date_list_ch0.append(date)
            ch0_num_apulse.append(ref_hist0.GetMean())
            ch0_num_apulse_err.append(ref_hist0.GetStdDev())
            #print(temp_ch0 / hist_ch0.GetEntries(), date)
        else:
            pass

        # If the histograms are not empty append the arrays
        if hist_ch1.GetEntries() > 0:
            num_apulse_ch1_list.append(100 * temp_ch1 / hist_ch0.GetEntries())
            date_list_ch1.append(date)
            ch1_num_apulse.append(ref_hist1.GetMean())
            ch1_num_apulse_err.append(ref_hist1.GetStdDev())
            #print(temp_ch1 / hist_ch1.GetEntries(), date)
        else:
            pass

        del hist_ch0
        del hist_ch1
        del file

    ch0_num_apulse_array = np.array(ch0_num_apulse)
    ch0_num_apulse_err_array = np.array(ch0_num_apulse_err)
    ch1_num_apulse_array = np.array(ch1_num_apulse)
    ch1_num_apulse_err_array = np.array(ch1_num_apulse_err)

    date_array_ch0 = np.array(date_list_ch0)
    date_array_ch1 = np.array(date_list_ch1)
    x0 = process_date(date_array_ch0)
    x1 = process_date(date_array_ch1)
    num_apulse_ch0_array = np.array(num_apulse_ch0_list)
    num_apulse_ch1_array = np.array(num_apulse_ch1_list)

    #plt.errorbar(x0, ch0_num_apulse_array, yerr=ch0_num_apulse_err_array, color='b')
    plt.plot(process_date(date_array_ch0), ch0_num_apulse_array, '.', color='b')
    plt.title("Average apulse number per waveform Channel 0")
    plt.xlabel("Exposure /days")
    plt.ylabel("Apulse per waveform")
    #plt.xscale('log')
    plt.ylim(0.2, 1)
    plt.grid(True)
    plt.show()

    #plt.errorbar(x0, ch0_num_apulse_array, yerr=ch0_num_apulse_err_array, color='b')
    plt.plot(process_date(date_array_ch1), ch1_num_apulse_array, '.', color='b')
    plt.title("Average apulse number per waveform Channel 1")
    plt.xlabel("Exposure /days")
    plt.ylabel("Apulse per waveform")
    #plt.xscale('log')
    plt.ylim(0.2, 1)
    # plt.ylim(2, 15)
    plt.grid(True)
    plt.show()

    plt.plot(process_exposure(date_array_ch0), num_apulse_ch0_array, '.', color='b')
    plt.title("Percentage of total apulses in [1500:2600]ns Channel 0")
    plt.xlabel("Exposure time atm-days")
    plt.ylabel("Percentage of apulses in sample range /%")
    plt.ylim(10, 50)
    #plt.xscale('log')
    plt.grid(True)
    plt.show()

    plt.plot(process_date(date_array_ch1), num_apulse_ch1_array, '.', color='b')
    plt.title("Percentage of total apulses in [1500:2600]ns Channel 1")
    plt.xlabel("Exposure time /days")
    plt.ylabel("Percentage of apulses in sample range /%")
    plt.ylim(10, 50)
    #plt.xscale('log')
    plt.grid(True)
    plt.show()


if __name__ == '__main__':
    main()
