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

    '''shape_2D_hist = file.Get("PMT_Pulse_Shapes")
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

    print("Total number of events: ", total)'''

    hist1 = file.Get("rise_time_mean;1")
    hist2 = file.Get("rise_time_stdev;1")
    hist3 = file.Get("rise_time_resolution;1")
    hist4 = file.Get("fall_time_mean;1")
    hist5 = file.Get("fall_time_stdev;1")
    hist6 = file.Get("fall_time_resolution;1")
    hist7 = file.Get("ratio_mean;1")
    hist8 = file.Get("ratio_stdev;1")
    hist9 = file.Get("ratio_resolution;1")

    hist_list = [hist1,hist2,hist3,hist4,hist5,hist6,hist7,hist8,hist9]

    for index, hist in enumerate(hist_list):
        name = hist.GetTitle()
        new_hist = ROOT.TH2F(name,name,
                             hist.GetNbinsX(),
                             0,
                             hist.GetNbinsX(),
                             hist.GetNbinsY(),
                             0,
                             hist.GetNbinsY())
        for i_x in range(0, hist.GetNbinsX() + 1):
            for i_y in range(0, hist.GetNbinsY() + 1):
                new_hist.Fill(i_x-1,i_y-1,round(hist.GetBinContent(i_x,i_y),1))
        c1 = ROOT.TCanvas(name)
        c1.cd()
        c1.SetGrid()
        new_hist.Draw("colztext")
        ROOT.gStyle.SetOptStat(0)
        c1.SaveAs("~/Desktop/PDFs/"+name+".pdf")
        del new_hist
        del c1

if __name__ == '__main__':
    main()
