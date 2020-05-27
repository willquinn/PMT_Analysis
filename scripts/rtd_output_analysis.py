import sys
sys.path.insert(1, '..')

import ROOT
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from functions.other_functions import io_parse_arguments, gaussian, chi2


def main():

    args = io_parse_arguments()
    input_data_filename = args.i
    output_data_filename = args.o

    run_number = input_data_filename.split("_")[1]

    try:
        data_file = open(input_data_filename, 'r')
    except FileNotFoundError as fnf_error:
        print(fnf_error)
        raise Exception("Error opening data file")

    # Create 2D histograms

    file = ROOT.TFile(input_data_filename, "READ")
    output_file = ROOT.TFile(output_data_filename, "RECREATE")
    file.cd()

    tree = file.T
    tree.Print()

    #for some reason:
    # row = column
    # column = ID number
    # OM_ID = row
    tree.Print()
    rise_time_hists = []
    fall_time_hists = []
    ratio_hists = []
    charge_hists = []

    ratio_mean_2D = ROOT.TH2F("ratio_mean_run_" + run_number, "ratio_mean_run_" + run_number,
                         20,0,20,13,0,13)
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
        temp1 = ROOT.TH1F("rise_OM_ID: " + str(i), "rise_OM_ID: "+ str(i), 100, 1800, 3500)
        temp2 = ROOT.TH1F("fall_OM_ID: " + str(i), "fall_OM_ID: " + str(i), 100, 1500, 1900)
        temp3 = ROOT.TH1F("ratio_OM_ID: " + str(i), "ratio_OM_ID: " + str(i), 100, 0.5, 2)
        temp4 = ROOT.TH1I("charge_OM_ID: " + str(i), "charge_OM_ID: " + str(i), 200, -100000, 0)
        #temp5 = ROOT.TH1I("amplitude_OM_ID: " + str(i), "charge_OM_ID: " + str(i), 2000, -2000, 0)
        rise_time_hists.append(temp1)
        fall_time_hists.append(temp2)
        ratio_hists.append(temp3)
        charge_hists.append(temp4)

    for event in tree:
        if event.column == 19:
            print("column 19")
        #print(event.event_num,event.OM_ID,event.row,event.column,event.amplitude,event.charge,event.baseline,event.rise_time,event.fall_time,event.peak_time)
        rise_time_hists[event.OM_ID].Fill(event.rise_time)
        fall_time_hists[event.OM_ID].Fill(event.fall_time)
        if event.fall_time == 0:
            pass
        else:
            ratio_hists[event.OM_ID].Fill(event.rise_time/event.fall_time)
        charge_hists[event.OM_ID].Fill(event.charge)

    output_file.cd()

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
    '''directory = file.GetDirectory("RTD calo histograms")
    directory.cd()

    hist_2d_mapping = ROOT.TH2I("Mapping_" + run_number,
                                "Mapping_" + run_number,
                                20, 0, 20,
                                13, 0, 13)

    hist_led_mu = ROOT.TH2D("mu_" + run_number,
                            "mu_" + run_number,
                            20, 0, 20,
                            13, 0, 13)

    hist_led_sig = ROOT.TH2D("sig_" + run_number,
                             "sig_" + run_number,
                             20, 0, 20,
                             13, 0, 13)

    hist_led_chi2 = ROOT.TH2D("chi2_" + run_number,
                              "chi2_" + run_number,
                              20, 0, 20,
                              13, 0, 13)

    for i_col in range(20):
        for i_row in range(13):
            print("Looking for: ", "hcharge_0." + str(i_col) + "." + str(i_row) + ";1")
            temp_hist = directory.Get("hcharge_0." + str(i_col) + "." + str(i_row) + ";1")

            try:
                num_entries = temp_hist.GetEntries()
            except:
                continue

            hist_2d_mapping.Fill(i_col, i_row, num_entries)

            temp_y_list = []
            temp_x_array = np.linspace(temp_hist.GetBinLowEdge(1), temp_hist.GetBinLowEdge(1) + int(float('{:.2f}'.format(temp_hist.GetBinWidth(2)*temp_hist.GetNbinsX()))), temp_hist.GetNbinsX())
            for i_bin in range(temp_hist.GetNbinsX()):
                temp_y_list.append(temp_hist.GetBinContent(i_bin))
            temp_y_array = np.array(temp_y_list)

            maximum = temp_hist.GetMaximum()

            peaks, _ = find_peaks(temp_y_array[1:], height=maximum/0.9, distance=75)

            if len(peaks) == 0 or len(peaks) > 2:
                color = 'tab:red'
                plt.plot(temp_x_array[1:], temp_y_array[1:], ".", color=color)
                plt.xlabel('Charge')
                plt.ylabel('Counts')
                plt.plot(temp_x_array[1:][peaks], temp_y_array[1:][peaks], "x", color='tab:green')
                plt.show(block=True)

            gaus_x_list = []
            gaus_y_list = []
            gaus_y_err_list = []

            for i_x in range(peaks[0] - 40, peaks[0] + 40):
                gaus_x_list.append(temp_x_array[1:][i_x])
                gaus_y_list.append(temp_y_array[1:][i_x])
                if temp_y_array[1:][i_x] == 0:
                    gaus_y_err_list.append(1.0)
                else:
                    gaus_y_err_list.append(np.sqrt(temp_y_array[1:][i_x]))

            gaus_x_array = np.array(gaus_x_list)
            gaus_y_array = np.array(gaus_y_list)

            x_array = np.linspace(temp_hist.GetBinLowEdge(1), temp_hist.GetBinLowEdge(1) + int(float('{:.2f}'.format(temp_hist.GetBinWidth(2) * temp_hist.GetNbinsX()))), 10000)

            p_guess = [-1, 1, 1, 1] # mu, sig, A, h
            p_bounds = [[-12, 0, 0, 0], [0, 10, 10000, 1]]

            popt, pcov = curve_fit(gaussian, gaus_x_list, gaus_y_list, p0=p_guess, bounds=p_bounds)

            plt.plot(temp_x_array[1:], temp_y_array[1:], ".", color='tab:blue')
            plt.plot(gaus_x_array, gaus_y_array, "x", color='tab:red')
            plt.xlabel('Charge')
            plt.ylabel('Counts')
            plt.plot(x_array[1:], gaussian(x_array[1:], *popt), "-", color='tab:green')
            plt.show(block=False)
            plt.pause(0.2)
            plt.close()

            hist_led_mu.Fill(i_col, i_row, float('{:.2f}'.format(popt[0])))
            hist_led_sig.Fill(i_col, i_row, float('{:.2f}'.format(popt[1])))

            #print(((np.array(gaus_y_list) - np.array(gaussian(gaus_x_list, *popt)))/np.array(gaus_y_err_list))**2)
            #print(np.sum(((np.array(gaus_y_list) - np.array(gaussian(gaus_x_list, *popt)))/np.array(gaus_y_err_list))**2)/(len(gaus_y_list) - len(popt)))

            cal_chi2 = chi2(gaus_y_list, gaus_y_err_list, gaussian(gaus_x_list, *popt), len(popt))

            hist_led_chi2.Fill(i_col, i_row, float('{:.2f}'.format(cal_chi2)))

        break

    c1 = ROOT.TCanvas("Mapping_" + run_number)
    c1.cd()
    c1.SetGrid()
    #new_shape_2D_hist.SetMinimum(0.95)
    #new_shape_2D_hist.SetMaximum(1)
    hist_2d_mapping.Draw("colztext")
    ROOT.gStyle.SetOptStat(0)
    c1.SaveAs("~/Desktop/Mapping" + run_number + ".pdf")
    del c1

    c2 = ROOT.TCanvas("led_mu_" + run_number)
    c2.cd()
    hist_led_mu.Draw("colztext")
    ROOT.gStyle.SetOptStat(0)

    c2.SaveAs("~/Desktop/led_mu_" + run_number + ".pdf")
    del c2

    c3 = ROOT.TCanvas("led_sig_" + run_number)
    c3.cd()
    hist_led_sig.Draw("colztext")
    ROOT.gStyle.SetOptStat(0)

    c3.SaveAs("~/Desktop/led_sig_" + run_number + ".pdf")
    del c3

    c4 = ROOT.TCanvas("led_chi2_" + run_number)
    c4.cd()
    hist_led_chi2.Draw("colztext")
    ROOT.gStyle.SetOptStat(0)

    c4.SaveAs("~/Desktop/led_chi2_" + run_number + ".pdf")
    del c4'''


if __name__ == '__main__':
    main()
