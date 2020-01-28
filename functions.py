import numpy as np


def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser(description="Input file names")
    parser.add_argument('-i', required=True, type=str, help='Input data file')
    parser.add_argument('-c', required=False, type=str, help='Config file')
    parser.add_argument('-sweep', required=False, type=str, help='Define whether you want to sweep. By default this is set to false')
    parser.add_argument('-r', required=False, type=str, help='Define whether you want ot recreate the file')
    args = parser.parse_args()
    return args


def get_normalisation_factor(vector: list):
    norm = 0.0
    for i in range(len(vector)):
        norm += vector[i] * vector[i]
    return np.sqrt(norm)


def get_date_time(input_data_file_name: str):
    temp1 = input_data_file_name.split("/")
    date = temp1[-2]
    temp2 = temp1[-1].split("_")
    temp3 = temp2[2].split("t")
    time = temp3[1]
    return date, time


def fit_bismuth_function_from_file(root_file_name: str):
    import ROOT
    root_file = ROOT.TFile(root_file_name, "UPDATE")

    directory_0 = root_file.GetDirectory("GAO607_000000_0000")
    directory_1 = root_file.GetDirectory("GAO612_000000_0000")

    spectra_0 = directory_0.Get("GAO607_000000_0000_pulse_charge_spectrum")
    spectra_1 = directory_1.Get("GAO612_000000_0000_pulse_charge_spectrum")

    spectra_list = []
    spectra_list.append(spectra_0)
    spectra_list.append(spectra_1)

    for i_pmt in range(len(spectra_list)):
        if i_pmt == 0:
            if spectra_list[i_pmt] is not 0x0:
                if spectra_list[i_pmt].GetEntries() == 0:
                    continue
                canvas_name = "GAO607_000000_0000"
                canvas = ROOT.TCanvas(canvas_name, canvas_name)
                fit = ROOT.TF1("fit",
                               "[0]*"
                               "(7.08*TMath::Gaus(x,[1],[2]) "
                               " + 1.84*TMath::Gaus(x,[1]*(1 + 72.144/975.651),[2]*1.036) "
                               " + 0.44*TMath::Gaus(x,[1]*(1 + 84.154/975.651),[2]*1.042)) "
                               " + 0.464*(exp(0.254*x)/(1 + exp((x - 28.43)/2.14)))",
                               33, 41)
                fit.SetParNames("A", "mu", "sigma")

                fit.SetParLimits(0, 0, 400)
                fit.SetParLimits(1, 34, 37)
                fit.SetParLimits(2, 0.8, 2)
                fit.SetParameters(319, 36, 1.09)

                spectra_list[i_pmt].Fit("fit", "", "", 33, 41)
                spectra_list[i_pmt].SetXTitle("Charge /pC")
                spectra_list[i_pmt].SetYTitle("Counts")
                spectra_list[i_pmt].SetTitle("Bi Integrated Charge Spectrum")

                resolution = (fit.GetParameter(2) / fit.GetParameter(1)) * 100.0
                chi2 = fit.GetChisquare() / fit.GetNDF()
                mu = fit.GetParameter(1)
                mu_err = fit.GetParError(1)
                sigma = fit.GetParameter(2)
                sigma_err = fit.GetParError(2)

                canvas.SetGrid()
                canvas.Update()

                canvas.Draw()
                ROOT.gStyle.SetOptStat(11111111)
                ROOT.gStyle.SetOptFit(1)
                ROOT.gStyle.SetStatY(0.9)
                ROOT.gStyle.SetStatX(0.9)
                ROOT.gStyle.SetStatW(0.8)
                ROOT.gStyle.SetStatH(0.1)
                canvas.SaveAs(canvas_name, "pdf")

        elif i_pmt == 1:
            if spectra_list[i_pmt] is not 0x0:
                if spectra_list[i_pmt].GetEntries() == 0:
                    continue
                canvas_name = "GAO612_000000_0000"
                canvas = ROOT.TCanvas(canvas_name, canvas_name)

                fit = ROOT.TF1("fit",
                               "[0]*"
                               "(7.08*TMath::Gaus(x,[1],[2]) "
                               " + 1.84*TMath::Gaus(x,[1]*(1 + 72.144/975.651),[2]*1.036) "
                               " + 0.44*TMath::Gaus(x,[1]*(1 + 84.154/975.651),[2]*1.042)) "
                               " + 0.515*(exp(0.2199*x)/(1 + exp((x - 31.68)/2.48)))",
                               37, 44)

                fit.SetParNames("A", "mu", "sigma")

                fit.SetParLimits(0, 0, 400)
                fit.SetParLimits(1, 37, 40)
                fit.SetParLimits(2, 0.8, 2)

                fit.SetParameters(319, 39, 1.09)

                spectra_list[i_pmt].Fit("fit", "", "", 37, 44)
                spectra_list[i_pmt].SetXTitle("Charge /pC")
                spectra_list[i_pmt].SetYTitle("Counts")
                spectra_list[i_pmt].SetTitle("Bi Integrated Charge Spectrum")

                resolution = (fit.GetParameter(2) / fit.GetParameter(1)) * 100.0
                chi2 = fit.GetChisquare() / fit.GetNDF()
                mu = fit.GetParameter(1)
                mu_err = fit.GetParError(1)
                sigma = fit.GetParameter(2)
                sigma_err = fit.GetParError(2)

                canvas.SetGrid()
                canvas.Update()

                canvas.Draw()
                ROOT.gStyle.SetOptStat(11111111)
                ROOT.gStyle.SetOptFit(1)
                ROOT.gStyle.SetStatY(0.9)
                ROOT.gStyle.SetStatX(0.9)
                ROOT.gStyle.SetStatW(0.8)
                ROOT.gStyle.SetStatH(0.1)
                canvas.SaveAs(canvas_name, "pdf")
