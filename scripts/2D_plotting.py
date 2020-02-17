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

    shape_2D_hist = file.Get("PMT_Pulse_Shapes")
    SD_2D_hist = file.Get("PMT_Pulse_Shapes_SD")
    map_hist = file.Get("PMT_Event_Mapping")

    new_shape_2D_hist = ROOT.TH2F("Pulse_Shapes",
                                  "Pulse_Shapes",
                                  shape_2D_hist.GetNbinsY(),
                                  0,
                                  shape_2D_hist.GetNbinsY(),
                                  shape_2D_hist.GetNbinsX(),
                                  0,
                                  shape_2D_hist.GetNbinsX())

    new_SD_2D_hist = ROOT.TH2F("Pulse_Shapes_SD",
                               "Pulse_Shapes_SD",
                               shape_2D_hist.GetNbinsY(),
                               0,
                               shape_2D_hist.GetNbinsY(),
                               shape_2D_hist.GetNbinsX(),
                               0,
                               shape_2D_hist.GetNbinsX())

    total = 0
    for i in range(0, shape_2D_hist.GetNbinsX()+1):
        for j in range(0, shape_2D_hist.GetNbinsY()+1):
            temp_shape = shape_2D_hist.GetBinContent(i, j)
            temp_SD = SD_2D_hist.GetBinContent(i, j)

            new_shape_2D_hist.SetBinContent(i, j, round(temp_shape, 3))
            new_SD_2D_hist.SetBinContent(i, j, round(temp_SD, 3))

            total += map_hist.GetBinContent(i, j)

    print("Total number of events: ", total)


    c1 = ROOT.TCanvas("PMT_Pulse_Shapes")
    c1.cd()
    c1.SetGrid()
    new_shape_2D_hist.SetMinimum(0.95)
    new_shape_2D_hist.SetMaximum(1)
    new_shape_2D_hist.Draw("colztext")
    ROOT.gStyle.SetOptStat(0)
    c1.SaveAs("~/Desktop/PMT_Pulse_Shapes.pdf")
    del c1

    c2 = ROOT.TCanvas("PMT_Pulse_Shapes_SD")
    c2.cd()
    new_SD_2D_hist.Draw("colztext")
    new_SD_2D_hist.SetMinimum(0)
    new_SD_2D_hist.SetMaximum(0.02)
    ROOT.gStyle.SetOptStat(0)

    c2.SaveAs("~/Desktop/PMT_Pulse_Shapes_SD.pdf")
    del c2

    c3 = ROOT.TCanvas("PMT_Event_Mapping")
    c3.cd()
    map_hist.Draw("colztext")
    ROOT.gStyle.SetOptStat(0)

    c3.SaveAs("~/Desktop/Map.pdf")
    del c3

if __name__ == '__main__':
    main()
