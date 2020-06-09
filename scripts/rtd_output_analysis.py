import sys
sys.path.insert(1, '..')

import ROOT
import tqdm
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from functions.other_functions import sncd_parse_arguments, gaussian, chi2
from scr.PMT_Array import PMT_Array
from scr.PMT_Waveform import PMT_Waveform


def create_templates():
    pass


def main():

    args = sncd_parse_arguments()
    input_data_filename = args.i
    output_data_filename = args.o
    config_file_name = args.c

    run_number = input_data_filename.split("_")[1]

    try:
        data_file = open(input_data_filename, 'r')
    except FileNotFoundError as fnf_error:
        print(fnf_error)
        raise Exception("Error opening data file")

    topology = [20, 13]
    num_pmts = topology[0]*topology[1]
    pmt_array = PMT_Array(topology, run_number)
    for i in range(topology[0]):
        for j in range(topology[1]):
            pmt_array.set_pmt_id("R:" + str(j) + "_C:" + str(i), j*20 + i)

    # Set the cuts you wish to apply
    # If you don't do this the defaults are used
    if config_file_name is not None:
        pmt_array.apply_setting(config_file_name)

    # Create 2D histograms

    file = ROOT.TFile(input_data_filename, "READ")
    file.cd()

    tree = file.T
    tree.Print()
    tree.Print()

    rise_times = [[] for i in range(num_pmts)]
    fall_times = [[] for i in range(num_pmts)]
    ratio_times = [[] for i in range(num_pmts)]
    charges = [[] for i in range(num_pmts)]
    amplitudes = [[] for i in range(num_pmts)]

    average_waveforms = [[0 for i in range(1064 - 200)] for i in range(num_pmts)]

    for event in tqdm.tqdm(tree):
        OM_ID = event.OM_ID
        rise_time = event.rise_time
        fall_time = event.fall_time
        charge = event.charge
        amplitude = event.amplitude
        waveform = event.waveform

        rise_times[OM_ID].append(rise_time)
        fall_times[OM_ID].append(fall_time)
        if fall_time == 0:
            pass
        else:
            ratio_times[OM_ID].append(rise_time/fall_time)
        charges[OM_ID].append(charge)
        amplitudes[OM_ID].append(amplitude)
        pmt_waveform = PMT_Waveform(waveform, pmt_array.get_pmt_object_number(OM_ID))

        for i in range(1064 - 200):
            # peak should be at 300
            average_waveforms[OM_ID][i] += pmt_waveform.get_pmt_waveform()[pmt_waveform.get_pmt_pulse_peak_position() - 100] - pmt_waveform.get_pmt_baseline()

        plt.plot(average_waveforms[OM_ID])
        plt.show()
        return

    output_file = ROOT.TFile(output_data_filename, "RECREATE")
    output_file.cd()

    rise_time_hists = []
    fall_time_hists = []
    ratio_hists = []
    charge_hists = []

    ratio_mean_2D = ROOT.TH2F("ratio_mean_run_" + run_number, "ratio_mean_run_" + run_number,
                              20, 0, 20, 13, 0, 13)
    ratio_stdev_2D = ROOT.TH2F("ratio_stdev_run_" + run_number, "ratio_stdev_run_" + run_number,
                               20, 0, 20, 13, 0, 13)
    rise_mean_2D = ROOT.TH2F("rise_mean_run_" + run_number, "rise_mean_run_" + run_number,
                             20, 0, 20, 13, 0, 13)
    rise_stdev_2D = ROOT.TH2F("rise_stdev_run_" + run_number, "rise_stdev_run_" + run_number,
                              20, 0, 20, 13, 0, 13)
    fall_mean_2D = ROOT.TH2F("fall_mean_run_" + run_number, "fall_mean_run_" + run_number,
                             20, 0, 20, 13, 0, 13)
    fall_stdev_2D = ROOT.TH2F("fall_stdev_run_" + run_number, "fall_stdev_run_" + run_number,
                              20, 0, 20, 13, 0, 13)
    for i in range(260):
        temp1 = ROOT.TH1F("rise_OM_ID: " + str(i), "rise_OM_ID: " + str(i), 100, 1800, 3500)
        temp2 = ROOT.TH1F("fall_OM_ID: " + str(i), "fall_OM_ID: " + str(i), 100, 1500, 1900)
        temp3 = ROOT.TH1F("ratio_OM_ID: " + str(i), "ratio_OM_ID: " + str(i), 100, 0.5, 2)
        temp4 = ROOT.TH1I("charge_OM_ID: " + str(i), "charge_OM_ID: " + str(i), 200, -10000, 0)
        # temp5 = ROOT.TH1I("amplitude_OM_ID: " + str(i), "charge_OM_ID: " + str(i), 2000, -2000, 0)
        rise_time_hists.append(temp1)
        fall_time_hists.append(temp2)
        ratio_hists.append(temp3)
        charge_hists.append(temp4)

    for i in range(260):
        row = int(i/20)
        column = i % 20
        ratio_mean_2D.Fill(column, row, ratio_hists[i].GetMean())
        ratio_stdev_2D.Fill(column, row, ratio_hists[i].GetStdDev())
        rise_mean_2D.Fill(column, row, rise_time_hists[i].GetMean())
        rise_stdev_2D.Fill(column, row, rise_time_hists[i].GetStdDev())
        fall_mean_2D.Fill(column, row, fall_time_hists[i].GetMean())
        fall_stdev_2D.Fill(column, row, fall_time_hists[i].GetStdDev())
        rise_time_hists[i].Write()
        fall_time_hists[i].Write()
        ratio_hists[i].Write()
        charge_hists[i].Write()

    ratio_mean_2D.Write()
    ratio_stdev_2D.Write()
    rise_mean_2D.Write()
    rise_stdev_2D.Write()
    fall_mean_2D.Write()
    fall_stdev_2D.Write()

    output_file.Close()


if __name__ == '__main__':
    main()
