import sys
sys.path.insert(1, '..')

import ROOT
import tqdm
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from functions.other_functions import sncd_parse_arguments, gaussian, chi2, inner_product
from scr.PMT_Array import PMT_Array
from scr.PMT_Waveform import PMT_Waveform


def get_rise_time(waveform, ratio):
    # waveform need baseline of zero
    peak = np.argmin(waveform)
    point = 0
    for i in range(peak, waveform.size):
        if waveform[i] > np.amin(waveform)*ratio:
            point = i
            break

    return point - peak


def create_templates(pmt_array: PMT_Array, tree, amp_cuts):
    average_waveforms = [[0 for i in range(690)] for j in range(pmt_array.get_pmt_total_number())]

    for event in tqdm.tqdm(tree):
        OM_ID = event.OM_ID
        waveform = event.waveform
        amplitude = (-1)*event.amplitude

        if amp_cuts[0][OM_ID] < amplitude < amp_cuts[1][OM_ID]:

            pmt_waveform = PMT_Waveform(waveform, pmt_array.get_pmt_object_number(OM_ID))

            temp = pmt_waveform.get_pmt_waveform()[pmt_waveform.get_pmt_pulse_peak_position() - 100:] - pmt_waveform.get_pmt_baseline()

            for i in range(len(average_waveforms[OM_ID])):
                average_waveforms[OM_ID][i] += temp[i]

    for om in range(len(average_waveforms)):
        norm = inner_product(np.array(average_waveforms[om]), np.array(average_waveforms[om]))
        if norm <= 0:
            average_waveforms[om] = np.array(average_waveforms[om])
        else:
            average_waveforms[om] = np.array(average_waveforms[om])/norm

    return average_waveforms


def write_average_waveforms(average_waveforms, write_file: str):

    with open(write_file, 'w') as f:
        for j in range(len(average_waveforms[0])):
            string = ''
            for i in range(len(average_waveforms)):
                string += str(average_waveforms[i][j]) + ','
            string += '\n'
            f.write(string)
    return


def read_average_waveforms(temp_file: str):
    average_waveforms = np.loadtxt(temp_file, delimiter=',', dtype={
        'names': ['{}'.format(i) for i in range(260)],
        'formats': ['f4' for i in range(260)]}, unpack=True)
    return average_waveforms


def main():

    # Handle the file inputs
    args = sncd_parse_arguments()
    input_data_filename = args.i
    output_data_filename = args.o
    config_file_name = args.c
    main_wall = args.w

    if main_wall == 'fr':
        main_wall = 1
    elif main_wall == 'it':
        main_wall = 0

    '''average_waveforms = read_average_waveforms("/Users/willquinn/Desktop/templates.txt")
    hist2D = ROOT.TH2I("", "",
              20, 0, 20, 13, 0, 13)
    for om in range(len(average_waveforms)):
        row = int(om / 20)
        column = om % 20
        if np.amax(average_waveforms[om]) > 0 or np.amin(average_waveforms[om]) < 0:
            hist2D.Fill(column,row,1)
        else:
            hist2D.Fill(column, row, 0)
    c1 = ROOT.TCanvas("")
    c1.cd()
    c1.SetGrid()
    hist2D.Draw("colz")
    ROOT.gStyle.SetOptStat(0)
    c1.SaveAs("~/Desktop/" + "template_test" + ".pdf")
    del c1

    return'''

    num_events = 1000

    # isolate the run number for naming convienience
    run_number = input_data_filename.split("_")[1]

    # Check to see if file exists - standard
    try:
        data_file = open(input_data_filename, 'r')
    except FileNotFoundError as fnf_error:
        print(fnf_error)
        raise Exception("Error opening data file")

    # Set up the pmt array
    topology = [20, 13]
    num_pmts = topology[0]*topology[1]
    pmt_array = PMT_Array(topology, run_number)
    for i in range(topology[0]):
        for j in range(topology[1]):
            num = i + topology[0]*j
            pmt_array.set_pmt_id("M:{}.{}.{}".format(main_wall, j, i), num)

    # Configure the array of PMTs - not as important here just yet
    if config_file_name is not None:
        pmt_array.apply_setting(config_file_name)

    # Open file
    file = ROOT.TFile(input_data_filename, "READ")
    file.cd()

    # The tree inside is called "T"
    tree = file.T

    # Create a series of lists to hold the information for manipulation
    rise_times = [[] for i in range(num_pmts)]
    fall_times = [[] for i in range(num_pmts)]
    new_rise_times = [[] for i in range(num_pmts)]
    ratio_times = [[] for i in range(num_pmts)]
    charges = [[] for i in range(num_pmts)]
    amplitudes = [[] for i in range(num_pmts)]

    uncut_amplitudes = [[] for i in range(num_pmts)]
    uncut_rise_times = [[] for i in range(num_pmts)]
    uncut_fall_times = [[] for i in range(num_pmts)]

    # Counter for testing
    event_counter = [0 for i in range(num_pmts)]

    # Run over file to create containers
    for event in tqdm.tqdm(tree):
        OM_ID = event.OM_ID
        charge = -1*event.charge
        amplitude = -1*event.amplitude/8
        fall_time = event.fall_time/256
        rise_time = event.rise_time/256

        charges[OM_ID].append(charge)
        amplitudes[OM_ID].append(amplitude)

        event_counter[OM_ID] += 1

        uncut_amplitudes[OM_ID].append(amplitude)
        uncut_fall_times[OM_ID].append(fall_time)
        uncut_rise_times[OM_ID].append(rise_time)

        if event_counter[0] == num_events:
            break

    # Handle the amplitudes to get a cut
    amp_hists = []
    max_amp = []
    amp_bin_edges = []
    bin_width = 30
    num_bins = 100
    my_bins = [i*bin_width for i in range(num_bins)]

    for ipmt in tqdm.tqdm(range(len(amplitudes))):
        hist,edges = np.histogram(amplitudes[ipmt], bins=my_bins)
        amp_hists.append(hist)
        max_amp.append(np.argmax(hist)*bin_width)
        amp_bin_edges.append(edges)

        '''plt.hist2d(uncut_amplitudes[ipmt],uncut_rise_times[ipmt], bins=(50,50), density=True, vmin=0, vmax=0.05)
        plt.xlabel("amplitude")
        plt.ylabel("rise time")
        plt.title(pmt_array.get_pmt_object_number(ipmt).get_pmt_id() + " Events:" + str(len(uncut_amplitudes[ipmt])))
        plt.colorbar()
        plt.savefig("/Users/willquinn/Desktop/PDFs/rise_vs_amp_om_" + str(ipmt))
        plt.close()
        plt.hist2d(uncut_amplitudes[ipmt], uncut_fall_times[ipmt], bins=(50,50), density=True, vmin=0, vmax=0.05)
        plt.xlabel("amplitude")
        plt.ylabel("fall time")
        plt.title(pmt_array.get_pmt_object_number(ipmt).get_pmt_id() + " Events:" + str(len(uncut_amplitudes[ipmt])))
        plt.colorbar()
        plt.savefig("/Users/willquinn/Desktop/PDFs/fall_vs_amp_om_" + str(ipmt))
        plt.close()
        plt.hist(uncut_amplitudes[ipmt], bins=num_bins, color="g", alpha=0.7, density=True)
        plt.xlabel("amplitude /mV")
        plt.ylabel("counts")
        plt.title(pmt_array.get_pmt_object_number(ipmt).get_pmt_id() + " Events:" + str(len(uncut_amplitudes[ipmt])))
        plt.savefig("/Users/willquinn/Desktop/PDFs/amp_om_" + str(ipmt))
        plt.close()'''


    # Counter for testing
    event_counter = [0 for i in range(num_pmts)]

    #average_waveforms = create_templates(pmt_array, tree, [np.array(max_amp) - bin_width*3, np.array(max_amp) + bin_width*4])
    #write_average_waveforms(average_waveforms, "/Users/willquinn/Desktop/templates.txt")

    for event in tqdm.tqdm(tree):
        OM_ID = event.OM_ID
        rise_time = event.rise_time/256
        fall_time = event.fall_time/256
        charge = -1*event.charge
        amplitude = -1*event.amplitude/8

        if max_amp[OM_ID] - bin_width*3 < amplitude < max_amp[OM_ID] + bin_width*4:
            waveform = PMT_Waveform(event.waveform, pmt_array.get_pmt_object_number(OM_ID))
            if waveform.get_pmt_waveform_reduced().size > 0:
                new_rise_time = get_rise_time(waveform.get_pmt_waveform_reduced(), 0.01)
                new_rise_times[OM_ID].append(new_rise_time)

            rise_times[OM_ID].append(rise_time)
            fall_times[OM_ID].append(fall_time)

            if fall_time > 0:
                ratio_times[OM_ID].append(rise_time/fall_time)

        event_counter[OM_ID] += 1
        if event_counter[0] == num_events:
            break

    output_file = ROOT.TFile(output_data_filename, "RECREATE")
    output_file.cd()

    # Grab the maximum value in each amp histogram
    ratio_mean_2D = ROOT.TH2F("ratio_mean_run_" + run_number, "ratio_mean_run_" + run_number,
                              20, 0, 20, 13, 0, 13)
    ratio_stdev_2D = ROOT.TH2F("ratio_stdev_run_" + run_number, "ratio_stdev_run_" + run_number,
                               20, 0, 20, 13, 0, 13)
    rise_mean_2D = ROOT.TH2F("rise_mean_run_" + run_number, "rise_mean_run_" + run_number,
                             20, 0, 20, 13, 0, 13)
    rise_stdev_2D = ROOT.TH2F("rise_stdev_run_" + run_number, "rise_stdev_run_" + run_number,
                              20, 0, 20, 13, 0, 13)
    new_rise_mean_2D = ROOT.TH2F("new_rise_mean_run_" + run_number, "new_rise_mean_run_" + run_number,
                             20, 0, 20, 13, 0, 13)
    new_rise_stdev_2D = ROOT.TH2F("new_rise_stdev_run_" + run_number, "new_rise_stdev_run_" + run_number,
                              20, 0, 20, 13, 0, 13)
    fall_mean_2D = ROOT.TH2F("fall_mean_run_" + run_number, "fall_mean_run_" + run_number,
                             20, 0, 20, 13, 0, 13)
    fall_stdev_2D = ROOT.TH2F("fall_stdev_run_" + run_number, "fall_stdev_run_" + run_number,
                              20, 0, 20, 13, 0, 13)
    amp_peak_2D = ROOT.TH2F("amp_peak_run_" + run_number, "amp_peak_run_" + run_number,
                              20, 0, 20, 13, 0, 13)

    for i in tqdm.tqdm(range(num_pmts)):
        plt.hist(np.array(rise_times[i])/np.array(rise_times[i]).size, bins=num_bins, color='b',alpha=0.8)
        plt.title("Rise times OM_ID:{} Mean:{} Std:{} Entries:{}".format(i,round(np.average(rise_times[i]),2),round(np.std(rise_times[i]),2),len(rise_times[i])))
        plt.savefig("/Users/willquinn/Desktop/PDFs/rise_time_om_"+str(i))
        plt.close()
        plt.hist(np.array(new_rise_times[i]) / np.array(new_rise_times[i]).size, bins=num_bins, color='b', alpha=0.8)
        plt.title("New Rise times OM_ID:{} Mean:{} Std:{} Entries:{}".format(i, round(np.average(new_rise_times[i]), 2),
                                                                         round(np.std(new_rise_times[i]), 2),
                                                                         len(new_rise_times[i])))
        plt.savefig("/Users/willquinn/Desktop/PDFs/new_rise_time_om_" + str(i))
        plt.close()
        plt.hist(np.array(fall_times[i])/np.array(fall_times[i]).size, bins=num_bins, color='b', alpha=0.8)
        plt.title("Fall times OM_ID:{} Mean:{} Std:{} Entries:{}".format(i, round(np.average(fall_times[i])), round(np.std(fall_times[i]),2),len(fall_times[i])))
        plt.savefig("/Users/willquinn/Desktop/PDFs/fall_time_om_" + str(i))
        plt.close()
        plt.hist(np.array(ratio_times[i])/np.array(ratio_times[i]).size, bins=num_bins, color='b', alpha=0.8)
        plt.title("Ratio rise/fall times OM_ID:{} Mean:{} Std:{} Entries:{}".format(i, round(np.average(ratio_times[i])), round(np.std(ratio_times[i]),2),len(ratio_times[i])))
        plt.savefig("/Users/willquinn/Desktop/PDFs/ratio_time_om_" + str(i))
        plt.close()

        row = int(i / 20)
        column = i % 20
        ratio_mean_2D.Fill(column, row, round(np.average(ratio_times[i]),2))
        ratio_stdev_2D.Fill(column, row, round(np.std(ratio_times[i]),2))
        rise_mean_2D.Fill(column, row, round(np.average(rise_times[i]),2))
        rise_stdev_2D.Fill(column, row, round(np.std(rise_times[i]),2))
        new_rise_mean_2D.Fill(column, row, round(np.average(new_rise_times[i]), 2))
        new_rise_stdev_2D.Fill(column, row, round(np.std(new_rise_times[i]), 2))
        fall_mean_2D.Fill(column, row, round(np.average(fall_times[i]),2))
        fall_stdev_2D.Fill(column, row, round(np.std(fall_times[i]),2))
        amp_peak_2D.Fill(column, row, max_amp[i])

    hists = [ratio_mean_2D,ratio_stdev_2D,rise_mean_2D,rise_stdev_2D,fall_stdev_2D,fall_mean_2D,amp_peak_2D, new_rise_mean_2D, new_rise_stdev_2D]
    for i in range(len(hists)):
        name = hists[i].GetName()
        c1 = ROOT.TCanvas(name)
        c1.cd()
        c1.SetGrid()
        hists[i].Draw("colz")
        ROOT.gStyle.SetOptStat(0)
        c1.SaveAs("~/Desktop/PDFs/" + name + ".pdf")
        del c1

    ratio_mean_2D.Write()
    ratio_stdev_2D.Write()
    rise_mean_2D.Write()
    rise_stdev_2D.Write()
    new_rise_mean_2D.Write()
    new_rise_stdev_2D.Write()
    fall_mean_2D.Write()
    fall_stdev_2D.Write()
    amp_peak_2D.Write()

    output_file.Close()

    '''OM_NUM = 0
    #plt.plot(amp_bin_edges[0][:-1],amp_hists[0],"r.")
    plt.hist(amplitudes[OM_NUM], bins=my_bins, color="g", alpha=0.8)
    plt.axvline(max_amp[OM_NUM] + bin_width/2, color="b", ls='--')
    plt.axvline(max_amp[OM_NUM] + bin_width/2 - bin_width*3, color="r", ls='--')
    plt.axvline(max_amp[OM_NUM] + bin_width/2 + bin_width*3, color="r", ls='--')
    plt.xlabel("amplitude")
    plt.ylabel("counts")
    plt.title("OM_ID:{}".format(OM_NUM))
    plt.yscale('log')
    plt.text(1,100,"Total events: {}".format(len(amplitudes[OM_NUM])))
    plt.text(1, 70, "Peak: {}".format(max_amp[OM_NUM] + bin_width/2))
    plt.grid()
    plt.show()'''

    '''x = np.array([(i + 100)/2.56 for i in range(len(average_waveforms[0]))])
    #colour_list0 = ['C0-','C1-','C2-','C3-','C4-']
    colour_list0 = ['C0-', 'C1-']
    #colour_list1 = ['C5-', 'C6-', 'C7-', 'C8-', 'C9-']
    colour_list1 = ['C5-', 'C6-']

    for i in range(len(colour_list0)):
        j = i
        print("OM_ID:",i + 40,"num:",event_counter[40 + i])
        print("OM_ID:",60 + i,"num:",event_counter[60 + i])
        while inner_product(np.array(average_waveforms[40 + i]), np.array(average_waveforms[40 + i])) <= 0:
            i += 1
        while inner_product(np.array(average_waveforms[60 + j]), np.array(average_waveforms[60 + j])) <= 0:
            j += 1

        plt.plot(x,-1*np.array(average_waveforms[40 + i])/inner_product(np.array(average_waveforms[40 + i]),np.array(average_waveforms[40 + i])),colour_list0[0],label=pmt_array.get_pmt_object_number(i + 40).get_pmt_id(),alpha=0.3)
        plt.plot(x,-1*np.array(average_waveforms[60+j]) / inner_product(np.array(average_waveforms[60+i]),np.array(average_waveforms[60+i])), colour_list1[0],
                 label=pmt_array.get_pmt_object_number(60+i).get_pmt_id(), alpha=0.7)
    plt.legend(loc="upper right")
    plt.xlim(np.amin(x), np.amax(x))
    plt.xlabel("time /ns")
    plt.ylabel("ADC /a.u")
    plt.yscale('log')
    plt.grid()
    plt.show()'''

    '''output_file = ROOT.TFile(output_data_filename, "RECREATE")
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
        charge_hists.append(temp4)'''


if __name__ == '__main__':
    main()
