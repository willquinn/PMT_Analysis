import sys

sys.path.insert(1, '..')

# import ROOT and bash commands
import ROOT

# import python plotting and numpy modules
import matplotlib.pyplot as plt
import numpy as np
import random as rand

# import stats module
from scipy.optimize import curve_fit

# import custom made classes
from functions.other_functions import sncd_parse_arguments, get_date_time, get_voltage, fit, chi2
from scr.PMT_Array import PMT_Array
from scr.PMT_Waveform import PMT_Waveform


def main():

    # Handle the input arguments:
    ##############################
    args = sncd_parse_arguments()
    input_data_file_name = args.i
    config_file_name = args.c
    output_file_name = args.o
    sweep_bool = args.sweep
    ##############################

    # Do some string manipulation to get the date and time from the file name
    #################################################
    date, time = get_date_time(input_data_file_name)
    voltage = int(get_voltage(input_data_file_name))
    #################################################

    # Check to see if the data file exists
    try:
        print(">>> Reading data from file: {}".format(input_data_file_name))
        date_file = open(input_data_file_name, 'r')
    except FileNotFoundError as fnf_error:
        print(fnf_error)
        raise Exception("Error opening data file {}".format(input_data_file_name))

    file = ROOT.TFile(input_data_file_name, "READ")
    output_file = ROOT.TFile(output_file_name, "RECREATE")
    file.cd()

    tree = file.T
    tree.Print()

    topology = [4, 1]
    pmt_array = PMT_Array(topology, date + "_" + time)
    pmt_array.set_pmt_id("GAO607_" + date + "_" + time, 0)
    pmt_array.set_pmt_id("GAO612_" + date + "_" + time, 1)
    pmt_array.set_pmt_id("Injected_GAO607_" + date + "_" + time, 2)
    pmt_array.set_pmt_id("Injected_GAO612_" + date + "_" + time, 3)

    # Set the cuts you wish to apply
    # If you don't do this the defaults are used
    if config_file_name is not None:
        pmt_array.apply_setting(config_file_name)

    if sweep_bool == "True":
        pmt_array.set_sweep_bool(True)
        if voltage == 1000:
            pmt_array.set_pmt_templates(
                "~/Desktop/191008_A1000_B1000_templates.root",
                ["A1000_B1000_Ch0_Template", "A1000_B1000_Ch1_Template"])
        elif voltage == 1400:
            pmt_array.set_pmt_templates(
                "~/Desktop/190621_A1400_B1400_templates.root",
                ["A1400_B1400_Ch0_Template", "A1400_B1400_Ch1_Template"])

    t_injected = [[], []]
    a_injected = [[], []]

    for event in tree:
        OM_ID = event.OM_ID
        event_num = event.event_num
        pulse_time = event.pulse_time
        pulse_amplitude = event.pulse_amplitude
        pulse_charge = event.pulse_charge
        pulse_baseline = event.pulse_baseline
        waveform = np.array(event.waveform)

        template_pulse = pmt_array.get_pmt_object_number(OM_ID).get_template_pmt_pulse()

        temp_amplitude = np.amin(template_pulse)
        temp_amp_pos = np.argmin(template_pulse)

        random_amp = rand.randrange(0, 101, 2)
        random_time = rand.randrange(800, pmt_array.get_pmt_object_number(OM_ID).get_setting("waveform_length") - template_pulse.size())

        t_injected[OM_ID].append(random_time)
        a_injected[OM_ID].append(random_amp)

        factor = random_amp / temp_amplitude

        injected_data = np.zeros_like(waveform)

        j = 0
        for i in range(injected_data.size()):
            if i >= random_time:
                injected_data[i] = template_pulse(j) * factor
                j += 1
                if j == template_pulse.size():
                    break

        pmt_waveform = PMT_Waveform(waveform, pmt_array.get_pmt_object_number(OM_ID))
        pmt_waveform_injected = PMT_Waveform(injected_data, pmt_array.get_pmt_object_number(OM_ID + 2))

    output_file.Close()


if __name__ == '__main__':
    main()
