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


def get_voltage(input_data_path: str):
    i = input_data_path.split("/")
    ii = i[-1]
    iii = ii.split("_")
    iv = iii[0]
    v = iv.split("A")
    vi = v[1]
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
        if date == 200214:
            temp_date = 100
        if date == 200215:
            temp_date = 101
        if date == 200216:
            temp_date = 102
        if date == 200217:
            temp_date = 103
        if date == 200218:
            temp_date = 104
        if date == 200219:
            temp_date = 105
        if date == 200220:
            temp_date = 106
        if date == 200220:
            temp_date = 107
        if date == 200221:
            temp_date = 108
        if date == 200222:
            temp_date = 109
        if date == 200223:
            temp_date = 110
        if date == 200224:
            temp_date = 111
        if date == 200225:
            temp_date = 112
        if date == 200226:
            temp_date = 113
        if date == 200227:
            temp_date = 114
        if date == 200228:
            temp_date = 115
        if date == 200229:
            temp_date = 116
        if date == 200301:
            temp_date = 117
        if date == 200302:
            temp_date = 118
        if date == 200303:
            temp_date = 119
        if date == 200304:
            temp_date = 120
        if date == 200305:
            temp_date = 121
        if date == 200306:
            temp_date = 122
        if date == 200307:
            temp_date = 123
        if date == 200308:
            temp_date = 124
        if date == 200309:
            temp_date = 125
        if date == 200310:
            temp_date = 126
        if date == 200311:
            temp_date = 127
        if date == 200312:
            temp_date = 128
        if date == 200313:
            temp_date = 129
        if date == 200314:
            temp_date = 130
        if date == 200315:
            temp_date = 131
        if date == 200316:
            temp_date = 132
        if date == 200317:
            temp_date = 133
        if date == 200318:
            temp_date = 134
        if date == 200319:
            temp_date = 135
        if date == 200320:
            temp_date = 136
        if date == 200321:
            temp_date = 137
        if date == 200322:
            temp_date = 138
        if date == 200323:
            temp_date = 139
        if date == 200324:
            temp_date = 140
        if date == 200325:
            temp_date = 141
        if date == 200326:
            temp_date = 142
        if date == 200327:
            temp_date = 143
        if date == 200328:
            temp_date = 144
        if date == 200329:
            temp_date = 145
        if date == 200330:
            temp_date = 146
        if date == 200331:
            temp_date = 147
        if date == 200401:
            temp_date = 148
        if date == 200402:
            temp_date = 149
        if date == 200403:
            temp_date = 150
        if date == 200404:
            temp_date = 151
        if date == 200405:
            temp_date = 152
        if date == 200406:
            temp_date = 153
        if date == 200407:
            temp_date = 154
        if date == 200408:
            temp_date = 155
        if date == 200409:
            temp_date = 156
        if date == 200410:
            temp_date = 157
        if date == 200411:
            temp_date = 158
        if date == 200412:
            temp_date = 159
        if date == 200413:
            temp_date = 160
        if date == 200414:
            temp_date = 161
        if date == 200415:
            temp_date = 162
        if date == 200416:
            temp_date = 163
        if date == 200417:
            temp_date = 164
        if date == 200418:
            temp_date = 165
        if date == 200419:
            temp_date = 166
        if date == 200420:
            temp_date = 167
        if date == 200421:
            temp_date = 168
        if date == 200422:
            temp_date = 169
        if date == 200423:
            temp_date = 170
        if date == 200424:
            temp_date = 171
        if date == 200425:
            temp_date = 172
        if date == 200426:
            temp_date = 173
        if date == 200427:
            temp_date = 174
        if date == 200428:
            temp_date = 175
        if date == 200429:
            temp_date = 176
        if date == 200430:
            temp_date = 177

        output_list.append(temp_date)

    output_array = np.array(output_list)
    assert output_array.size == date_array.size
    return output_array


def process_exposure(date_array: np.array):
    output_list = []
    for index, date in enumerate(date_array):
        temp_date = 0
        if date == 191008:
            temp_date = 0/1000000
        if date == 191009:
            temp_date = 1/1000000
        if date == 191010:
            temp_date = 2/1000000
        if date == 191011:
            temp_date = 3/1000000
        if date == 191012:
            temp_date = 4/1000000
        if date == 191013:
            temp_date = 5/1000000
        if date == 191014:
            temp_date = 6/1000000
        if date == 191015:
            temp_date = 7/1000000
        if date == 191016:
            temp_date = 8/1000000
        if date == 191017:
            temp_date = 9/1000000
        if date == 191018:
            temp_date = 10/1000000
        if date == 191019:
            temp_date = 11/1000000
        if date == 191020:
            temp_date = 12/1000000
        if date == 191021:
            temp_date = 13/1000000
        if date == 191022:
            temp_date = 14/1000000
        if date == 191023:
            temp_date = 15/1000000
        if date == 191024:
            temp_date = 16/1000000
        if date == 191025:
            temp_date = 17/1000000
        if date == 191026:
            temp_date = 18/1000000
        if date == 191027:
            temp_date = 19/1000000
        if date == 191028:
            temp_date = 20/1000000
        if date == 191029:
            temp_date = 21/1000000
        if date == 191030:
            temp_date = 22/1000000
        if date == 191031:
            temp_date = 23/1000000
        if date == 191101:
            temp_date = 24/1000000
        if date == 191102:
            temp_date = 25/1000000
        if date == 191103:
            temp_date = 26/1000000
        if date == 191104:
            temp_date = 27/1000000
        if date == 191105:
            temp_date = 28/1000000
        if date == 191106:
            temp_date = 29/1000000
        if date == 191107:
            temp_date = 1/100
        if date == 191108:
            temp_date = 2/100
        if date == 191109:
            temp_date = 3/100
        if date == 191110:
            temp_date = 4/100
        if date == 191111:
            temp_date = 5/100
        if date == 191112:
            temp_date = 6/100
        if date == 191113:
            temp_date = 7/100
        if date == 191114:
            temp_date = 8/100
        if date == 191115:
            temp_date = 9/100
        if date == 191116:
            temp_date = 10/100
        if date == 191117:
            temp_date = 11/100
        if date == 191118:
            temp_date = 12/100
        if date == 191119:
            temp_date = 13/100
        if date == 191120:
            temp_date = 14/100
        if date == 191121:
            temp_date = 15/100
        if date == 191122:
            temp_date = 16/100
        if date == 191123:
            temp_date = 17/100
        if date == 191124:
            temp_date = 18/100
        if date == 191125:
            temp_date = 19/100
        if date == 191126:
            temp_date = 20/100
        if date == 191127:
            temp_date = 21/100
        if date == 191128:
            temp_date = 22/100
        if date == 191129:
            temp_date = 23/100
        if date == 191130:
            temp_date = 24/100
        if date == 191201:
            temp_date = 25/100
        if date == 191202:
            temp_date = 26/100
        if date == 191203:
            temp_date = 27/100
        if date == 191204:
            temp_date = 28/100
        if date == 191205:
            temp_date = 29/100
        if date == 191206:
            temp_date = 30/100
        if date == 191207:
            temp_date = 31/100
        if date == 191208:
            temp_date = 32/100
        if date == 191209:
            temp_date = 33/100
        if date == 191210:
            temp_date = 34/100
        if date == 191211:
            temp_date = 35/100
        if date == 191212:
            temp_date = 36/100
        if date == 191213:
            temp_date = 37/100
        if date == 191214:
            temp_date = 38/100
        if date == 191215:
            temp_date = 39/100
        if date == 191216:
            temp_date = 40/100
        if date == 191217:
            temp_date = 41/100
        if date == 191218:
            temp_date = 42/100
        if date == 191219:
            temp_date = 43/100
        if date == 191220:
            temp_date = 44/100
        if date == 191221:
            temp_date = 45/100
        if date == 191222:
            temp_date = 46/100
        if date == 191223:
            temp_date = 47/100
        if date == 191224:
            temp_date = 48/100
        if date == 191225:
            temp_date = 49/100
        if date == 191226:
            temp_date = 50/100
        if date == 191227:
            temp_date = 51/100
        if date == 191228:
            temp_date = 52/100
        if date == 191229:
            temp_date = 53/100
        if date == 191230:
            temp_date = 54/100
        if date == 191231:
            temp_date = 55/100
        if date == 200101:
            temp_date = 56/100
        if date == 200102:
            temp_date = 57/100
        if date == 200103:
            temp_date = 58/100
        if date == 200104:
            temp_date = 59/100
        if date == 200105:
            temp_date = 60/100
        if date == 200106:
            temp_date = 61/100
        if date == 200107:
            temp_date = 62/100
        if date == 200108:
            temp_date = 63/100
        if date == 200109:
            temp_date = 64/100
        if date == 200110:
            temp_date = 65/100
        if date == 200111:
            temp_date = 66/100
        if date == 200112:
            temp_date = 67/100
        if date == 200113:
            temp_date = 68/100
        if date == 200114:
            temp_date = 69/100
        if date == 200115:
            temp_date = 70/100
        if date == 200116:
            temp_date = 71/100
        if date == 200117:
            temp_date = 72/100
        if date == 200118:
            temp_date = 73/100
        if date == 200119:
            temp_date = 74/100
        if date == 200120:
            temp_date = 75/100
        if date == 200121:
            temp_date = 76/100
        if date == 200122:
            temp_date = 77/100
        if date == 200123:
            temp_date = 78/100
        if date == 200124:
            temp_date = 79/100
        if date == 200125:
            temp_date = 80/100
        if date == 200126:
            temp_date = 81/100
        if date == 200127:
            temp_date = 82/100
        if date == 200128:
            temp_date = 83/100
        if date == 200129:
            temp_date = 84/100
        if date == 200130:
            temp_date = 85/100
        if date == 200131:
            temp_date = 86/100
        if date == 200201:
            temp_date = 87/100
        if date == 200202:
            temp_date = 88/100
        if date == 200203:
            temp_date = 89/100
        if date == 200204:
            temp_date = 90/100
        if date == 200205:
            temp_date = 91/100
        if date == 200206:
            temp_date = 92/100
        if date == 200207:
            temp_date = 93/100
        if date == 200208:
            temp_date = 94/100
        if date == 200209:
            temp_date = 95/100
        if date == 200210:
            temp_date = 96/100
        if date == 200211:
            temp_date = 97/100
        if date == 200212:
            temp_date = 1/10 + 97/100
        if date == 200213:
            temp_date = 2/10 + 97/100
        if date == 200214:
            temp_date = 3/10 + 97/100
        if date == 200215:
            temp_date = 4/10 + 97/100
        if date == 200216:
            temp_date = 5/10 + 97/100
        if date == 200217:
            temp_date = 6/10 + 97/100
        if date == 200218:
            temp_date = 7/10 + 97/100
        if date == 200219:
            temp_date = 8/10 + 97/100
        if date == 200220:
            temp_date = 9/10 + 97/100
        if date == 200220:
            temp_date = 10/10 + 97/100
        if date == 200221:
            temp_date = 11/10 + 97/100
        if date == 200222:
            temp_date = 12/10 + 97/100
        if date == 200223:
            temp_date = 13/10 + 97/100
        if date == 200224:
            temp_date = 14/10 + 97/100
        if date == 200225:
            temp_date = 15/10 + 97/100
        if date == 200226:
            temp_date = 16/10 + 97/100
        if date == 200227:
            temp_date = 17/10 + 97/100
        if date == 200228:
            temp_date = 18/10 + 97/100
        if date == 200229:
            temp_date = 19/10 + 97/100
        if date == 200301:
            temp_date = 20/10 + 97/100
        if date == 200302:
            temp_date = 21/10 + 97/100
        if date == 200303:
            temp_date = 22/10 + 97/100
        if date == 200304:
            temp_date = 23/10 + 97/100
        if date == 200305:
            temp_date = 24/10 + 97/100
        if date == 200306:
            temp_date = 25/10 + 97/100
        if date == 200307:
            temp_date = 26/10 + 97/100
        if date == 200308:
            temp_date = 27/10 + 97/100
        if date == 200309:
            temp_date = 28/10 + 97/100
        if date == 200310:
            temp_date = 29/10 + 97/100
        if date == 200311:
            temp_date = 30/10 + 97/100
        if date == 200312:
            temp_date = 31/10 + 97/100
        if date == 200313:
            temp_date = 32/10 + 97/100
        if date == 200314:
            temp_date = 33/10 + 97/100
        if date == 200315:
            temp_date = 34/10 + 97/100
        if date == 200316:
            temp_date = 35/10 + 97/100
        if date == 200317:
            temp_date = 36/10 + 97/100
        if date == 200318:
            temp_date = 37/10 + 97/100
        if date == 200319:
            temp_date = 38/10 + 97/100
        if date == 200320:
            temp_date = 39/10 + 97/100
        if date == 200321:
            temp_date = 40/10 + 97/100
        if date == 200322:
            temp_date = 41/10 + 97/100
        if date == 200323:
            temp_date = 42/10 + 97/100
        if date == 200324:
            temp_date = 43/10 + 97/100
        if date == 200325:
            temp_date = 44/10 + 97/100
        if date == 200326:
            temp_date = 45/10 + 97/100
        if date == 200327:
            temp_date = 46/10 + 97/100
        if date == 200328:
            temp_date = 47/10 + 97/100
        if date == 200329:
            temp_date = 48/10 + 97/100
        if date == 200330:
            temp_date = 49/10 + 97/100
        if date == 200331:
            temp_date = 50/10 + 97/100
        if date == 200401:
            temp_date = 51/10 + 97/100
        if date == 200402:
            temp_date = 52/10 + 97/100
        if date == 200403:
            temp_date = 53/10 + 97/100
        if date == 200404:
            temp_date = 54/10 + 97/100
        if date == 200405:
            temp_date = 55/10 + 97/100
        if date == 200406:
            temp_date = 56/10 + 97/100
        if date == 200407:
            temp_date = 57/10 + 97/100
        if date == 200408:
            temp_date = 58/10 + 97/100
        if date == 200409:
            temp_date = 59/10 + 97/100
        if date == 200410:
            temp_date = 60/10 + 97/100
        if date == 200411:
            temp_date = 61/10 + 97/100
        if date == 200412:
            temp_date = 62/10 + 97/100
        if date == 200413:
            temp_date = 63/10 + 97/100
        if date == 200414:
            temp_date = 64/10 + 97/100
        if date == 200415:
            temp_date = 65/10 + 97/100
        if date == 200416:
            temp_date = 66/10 + 97/100
        if date == 200417:
            temp_date = 67/10 + 97/100
        if date == 200418:
            temp_date = 68/10 + 97/100
        if date == 200419:
            temp_date = 69/10 + 97/100
        if date == 200420:
            temp_date = 70/10 + 97/100
        if date == 200421:
            temp_date = 71/10 + 97/100
        if date == 200422:
            temp_date = 72/10 + 97/100
        if date == 200423:
            temp_date = 73/10 + 97/100
        if date == 200424:
            temp_date = 74/10 + 97/100
        if date == 200425:
            temp_date = 75/10 + 97/100
        if date == 200426:
            temp_date = 76/10 + 97/100
        if date == 200427:
            temp_date = 77/10 + 97/100
        if date == 200428:
            temp_date = 78/10 + 97/100
        if date == 200429:
            temp_date = 79/10 + 97/100
        if date == 200430:
            temp_date = 80/10 + 97/100

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


def linear(x, m, c):
    y = []
    for i in range(len(x)):
        y.append(m*x[i] + c)
    return y


def fit_0(x, mu, sig, a):
    y = []
    for i in range(len(x)):
        calc =  a * (
                ( 7.08*gaus(x[i], mu, sig ) +
                ( 1.84*gaus(x[i], mu*(1 + 72.144/975.651), sig*1.036) ) +
                ( 0.44*gaus(x[i], mu*(1 + 84.154/975.651), sig*1.042) ) ) +
                0.464 *
                ( np.exp( 0.254*x[i] )/( 1 + np.exp( ( x[i] - 28.43 )/2.14) ) )
        )
        y.append(calc)
    return y


def fit_1(x, mu, sig, a):
    y = []
    for i in range(len(x)):
        calc = a * (
                (7.08 * gaus(x[i], mu, sig) +
                 (1.84 * gaus(x[i], mu * (1 + 72.144 / 975.651), sig * 1.036)) +
                 (0.44 * gaus(x[i], mu * (1 + 84.154 / 975.651), sig * 1.042))) +
                0.515 *
                (np.exp(0.2199 * x[i]) / (1 + np.exp((x[i] - 31.68) / 2.48)))
        )
        y.append(calc)
    return y


def fit(x, mu, sig, a, b, c, d, e):
    y = []
    for i in range(len(x)):
        calc = a * (
                (7.08 * gaus(x[i], mu, sig) +
                 (1.84 * gaus(x[i], mu * (1 + 72.144 / 975.651), sig * 1.036)) +
                 (0.44 * gaus(x[i], mu * (1 + 84.154 / 975.651), sig * 1.042))) +
                b *
                (np.exp(c * x[i]) / (1 + np.exp((x[i] - d) / e)))
        )
        y.append(calc)
    return y


def chi2(y_obs, y_err, y_exp, n_par):
    chi2 = 0
    ndof = len(y_obs) - n_par - 1
    for i in range(len(y_exp)):
        chi2 += ((y_exp[i] - y_obs[i])/y_err[i])**2
    chi2 = chi2/ndof
    return chi2


def gaus(x, mu, sig):
    return np.exp(-( ( ( x - mu )/sig )**2 )/2)/( np.sqrt( 2*np.pi )*sig )

