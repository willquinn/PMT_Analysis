import sys
sys.path.insert(1, '..')

import ROOT
import numpy as np
from functions.other_functions import io_parse_arguments


def main():
    args = io_parse_arguments()
    input_data_filename = args.i

    try:
        data_file = open(input_data_filename, 'r')
    except FileNotFoundError as fnf_error:
        print(fnf_error)
        raise Exception("Error opening data file")

    file = ROOT.TFile(input_data_filename, "READ")
    baseline_2D_hist = ROOT.TH2F("baseline_average",
                                 "baseline_average",
                                 20,0,20,13,0,13)
    baseline_sd_2D_hist = ROOT.TH2F("baseline_sd",
                                 "baseline_sd",
                                 20, 0, 20, 13, 0, 13)

    for i_pmt in range(260):
        directory = file.GetDirectory("{}".format(i_pmt))
        if 0 <= i_pmt <= 19:
            i_row = 0
        elif 13 <= i_pmt <= 39:
            i_row = 1
        elif 13 <= i_pmt <= 59:
            i_row = 2
        elif 13 <= i_pmt <= 79:
            i_row = 3
        elif 13 <= i_pmt <= 99:
            i_row = 4
        elif 13 <= i_pmt <= 119:
            i_row = 5
        elif 13 <= i_pmt <= 139:
            i_row = 6
        elif 13 <= i_pmt <= 159:
            i_row = 7
        elif 13 <= i_pmt <= 179:
            i_row = 8
        elif 13 <= i_pmt <= 199:
            i_row = 9
        elif 13 <= i_pmt <= 219:
            i_row = 10
        elif 13 <= i_pmt <= 239:
            i_row = 11
        elif 13 <= i_pmt <= 259:
            i_row = 12
        i_col = i_pmt - i_row*20

        try:
            mean = round(directory.Get("{}_baseline".format(i_pmt)).GetMean(),1)
            sd = round(directory.Get("{}_baseline".format(i_pmt)).GetStdDev(),3)
        except:
            continue

        baseline_2D_hist.Fill(i_col, i_row, mean)
        baseline_sd_2D_hist.Fill(i_col, i_row, sd)

        del directory

    c1 = ROOT.TCanvas("baseline_average")
    c1.cd()
    c1.SetGrid()
    baseline_2D_hist.SetMinimum(2046)
    baseline_2D_hist.SetMaximum(2049)
    baseline_2D_hist.Draw("colztext")
    ROOT.gStyle.SetOptStat(0)
    c1.SaveAs("~/Desktop/baseline_average.pdf")
    del c1

    c2 = ROOT.TCanvas("baseline_sd")
    c2.cd()
    baseline_sd_2D_hist.Draw("colztext")
    #new_SD_2D_hist.SetMinimum(0)
    #new_SD_2D_hist.SetMaximum(0.02)
    ROOT.gStyle.SetOptStat(0)

    c2.SaveAs("~/Desktop/baseline_sd.pdf")
    del c2


if __name__ == '__main__':
    main()