import numpy as np


def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser(description="Input file names")
    parser.add_argument('-i', required=True, type=str, help='Input data file')
    parser.add_argument('-c', required=False, type=str, help='Config file')
    parser.add_argument('-sweep', required=False, type=str, help='Define whether you want to sweep. By default this is set to false')
    parser.add_argument('-r', required=False, type=str, help='Define whether you want ot recreate the file')
    parser.add_argument('-f', required=False, type=str, help='Define if you want the Bismuth spectrum to be created')
    parser.add_argument('-t', required=False, type=list, help='Define the topology of the PMT Array')
    args = parser.parse_args()
    return args


def sncd_parse_arguments():
    import argparse
    parser = argparse.ArgumentParser(description="Input file names")
    parser.add_argument('-i', required=True, type=str, help='Input data file path')
    parser.add_argument('-c', required=False, type=str, help='Config file')
    parser.add_argument('-sweep', required=False, type=str, help='Define whether you want to sweep. By default this is set to false')
    parser.add_argument('-topology', required=False, type=list, help='Define the topology of the PMT Array')
    parser.add_argument('-t', required=False, type=str, help='Path for the templates')
    parser.add_argument('-o', required=True, type=str, help='Output file name')
    args = parser.parse_args()
    return args


def io_parse_arguments():
    import argparse
    parser = argparse.ArgumentParser(description="Input file names")
    parser.add_argument('-i', required=True, type=str, help='Input data file path')
    parser.add_argument('-o', required=False, type=str, help='Output data file path')
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


def get_run_number(input_data_path: str):
    i = input_data_path.split("/")
    ii = i[-1]
    iii = ii.split("_")
    iv = iii[2]
    v = iv.split(".")
    vi = v[0]
    return str(vi.strip())


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
                canvas_name = "GAO612_000000_0000.pdf"
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


def get_data_path(input_data_file: str):
    return input_data_file.split('filenames')[0]

def process_date(date_array: np.array):
    output_list = []
    for index, date in enumerate(date_array):
        temp_date = 0
        if date == 191008:
            temp_date = -29
        if date == 191009:
            temp_date = -28
        if date == 191010:
            temp_date = -27
        if date == 191011:
            temp_date = -26
        if date == 191012:
            temp_date = -25
        if date == 191013:
            temp_date = -24
        if date == 191014:
            temp_date = -23
        if date == 191015:
            temp_date = -22
        if date == 191016:
            temp_date = -21
        if date == 191017:
            temp_date = -20
        if date == 191018:
            temp_date = -19
        if date == 191019:
            temp_date = -18
        if date == 191020:
            temp_date = -17
        if date == 191021:
            temp_date = -16
        if date == 191022:
            temp_date = -15
        if date == 191023:
            temp_date = -14
        if date == 191024:
            temp_date = -13
        if date == 191025:
            temp_date = -12
        if date == 191026:
            temp_date = -11
        if date == 191027:
            temp_date = -10
        if date == 191028:
            temp_date = -9
        if date == 191029:
            temp_date = -8
        if date == 191030:
            temp_date = -7
        if date == 191031:
            temp_date = -6
        if date == 191101:
            temp_date = -5
        if date == 191102:
            temp_date = -4
        if date == 191103:
            temp_date = -3
        if date == 191104:
            temp_date = -2
        if date == 191105:
            temp_date = -1
        if date == 191106:
            temp_date = 0
        if date == 191107:
            temp_date = 1
        if date == 191108:
            temp_date = 2
        if date == 191109:
            temp_date = 3
        if date == 191110:
            temp_date = 4
        if date == 191111:
            temp_date = 5
        if date == 191112:
            temp_date = 6
        if date == 191113:
            temp_date = 7
        if date == 191114:
            temp_date = 8
        if date == 191115:
            temp_date = 9
        if date == 191116:
            temp_date = 10
        if date == 191117:
            temp_date = 11
        if date == 191118:
            temp_date = 12
        if date == 191119:
            temp_date = 13
        if date == 191120:
            temp_date = 14
        if date == 191121:
            temp_date = 15
        if date == 191122:
            temp_date = 16
        if date == 191123:
            temp_date = 17
        if date == 191124:
            temp_date = 18
        if date == 191125:
            temp_date = 19
        if date == 191126:
            temp_date = 20
        if date == 191127:
            temp_date = 21
        if date == 191128:
            temp_date = 22
        if date == 191129:
            temp_date = 23
        if date == 191130:
            temp_date = 24
        if date == 191201:
            temp_date = 25
        if date == 191202:
            temp_date = 26
        if date == 191203:
            temp_date = 27
        if date == 191204:
            temp_date = 28
        if date == 191205:
            temp_date = 29
        if date == 191206:
            temp_date = 30
        if date == 191207:
            temp_date = 31
        if date == 191208:
            temp_date = 32
        if date == 191209:
            temp_date = 33
        if date == 191210:
            temp_date = 34
        if date == 191211:
            temp_date = 35
        if date == 191212:
            temp_date = 36
        if date == 191213:
            temp_date = 37
        if date == 191214:
            temp_date = 38
        if date == 191215:
            temp_date = 39
        if date == 191216:
            temp_date = 40
        if date == 191217:
            temp_date = 41
        if date == 191218:
            temp_date = 42
        if date == 191219:
            temp_date = 43
        if date == 191220:
            temp_date = 44
        if date == 191221:
            temp_date = 45
        if date == 191222:
            temp_date = 46
        if date == 191223:
            temp_date = 47
        if date == 191224:
            temp_date = 48
        if date == 191225:
            temp_date = 49
        if date == 191226:
            temp_date = 50
        if date == 191227:
            temp_date = 51
        if date == 191228:
            temp_date = 52
        if date == 191229:
            temp_date = 53
        if date == 191230:
            temp_date = 54
        if date == 191231:
            temp_date = 55
        if date == 200101:
            temp_date = 56
        if date == 200102:
            temp_date = 57
        if date == 200103:
            temp_date = 58
        if date == 200104:
            temp_date = 59
        if date == 200105:
            temp_date = 60
        if date == 200106:
            temp_date = 61
        if date == 200107:
            temp_date = 62
        if date == 200108:
            temp_date = 63
        if date == 200109:
            temp_date = 64
        if date == 200110:
            temp_date = 65
        if date == 200111:
            temp_date = 66
        if date == 200112:
            temp_date = 67
        if date == 200113:
            temp_date = 68
        if date == 200114:
            temp_date = 69
        if date == 200115:
            temp_date = 70
        if date == 200116:
            temp_date = 71
        if date == 200117:
            temp_date = 72
        if date == 200118:
            temp_date = 73
        if date == 200119:
            temp_date = 74
        if date == 200120:
            temp_date = 75
        if date == 200121:
            temp_date = 76
        if date == 200122:
            temp_date = 77
        if date == 200123:
            temp_date = 78
        if date == 200124:
            temp_date = 79
        if date == 200125:
            temp_date = 80
        if date == 200126:
            temp_date = 81
        if date == 200127:
            temp_date = 82
        if date == 200128:
            temp_date = 83
        if date == 200129:
            temp_date = 84
        if date == 200130:
            temp_date = 85
        if date == 200131:
            temp_date = 86
        if date == 200201:
            temp_date = 87
        if date == 200202:
            temp_date = 88
        if date == 200203:
            temp_date = 89
        if date == 200204:
            temp_date = 90
        if date == 200205:
            temp_date = 91
        if date == 200206:
            temp_date = 92
        if date == 200207:
            temp_date = 93
        if date == 200208:
            temp_date = 94
        if date == 200209:
            temp_date = 95
        if date == 200210:
            temp_date = 96
        if date == 200211:
            temp_date = 97
        if date == 200212:
            temp_date = 98
        if date == 200213:
            temp_date = 99

        output_list.append(temp_date)

    output_array = np.array(output_list)
    assert output_array.size == date_array.size
    return output_array


def chi2(y_obs, y_err, y_exp, n_par):
    chi2 = 0
    ndof = len(y_obs) - n_par - 1
    for i in range(len(y_exp)):
        chi2 += ((y_exp[i] - y_obs[i])/y_err[i])**2
    chi2 = chi2/ndof
    return chi2


def gaussian(x, mean, sigma, amplitude, height):
    y = []
    for i in range(len(x)):
        y.append(amplitude*np.exp((-(x[i] - mean)**2)/(2*sigma**2)) + height)
    return y
