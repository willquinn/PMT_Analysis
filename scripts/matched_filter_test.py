import sys

sys.path.insert(1, '..')

# import ROOT and bash commands
import ROOT
import progressbar
import tqdm

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
    sweep_bool = "True"
    ##############################

    # Do some string manipulation to get the date and time from the file name
    #################################################
    date, time = "0000", "0000"
    voltage = 1400
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
                ["A1000_B1000_Ch0_Template", "A1000_B1000_Ch1_Template", "A1000_B1000_Ch0_Template", "A1400_B1400_Ch1_Template"])
        elif voltage == 1400:
            pmt_array.set_pmt_templates(
                "~/Desktop/190621_A1400_B1400_templates.root",
                ["A1400_B1400_Ch0_Template", "A1400_B1400_Ch1_Template", "A1400_B1400_Ch0_Template", "A1400_B1400_Ch1_Template"])

    t_injected = [[], []]
    a_injected = [[], []]

    num = [0, 0]
    enum = [0, 0]
    iterator = 0

    injected_array = [[[], []], [[], []]]
    best_array = [[[], []], [[], []]]
    random_array = [[[], []], [[], []]]
    x_best_array = [[], []]

    #bar = progressbar.ProgressBar(maxval=tree.GetEntries(), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    #bar.start()

    for event in tree:
        OM_ID = event.OM_ID
        waveform = np.array(event.waveform)

        # Get template for injecting into data
        template_pulse = pmt_array.get_pmt_object_number(OM_ID).get_template_pmt_pulse()

        # Get the amplitude of the template so we can normalise to 1
        temp_amplitude = np.amin(template_pulse)

        # Create random numbers
        random_amp = rand.randrange(0, 101, 2)
        random_time = rand.randrange(800, pmt_array.get_pmt_object_number(OM_ID).get_waveform_length() - template_pulse.size)

        '''print("A: ",random_amp)
        print("t: ",random_time)'''

        factor = random_amp / temp_amplitude

        injected_data = np.zeros_like(waveform)

        j = 0
        for i in range(injected_data.size):
            if i >= random_time:
                injected_data[i] = int(template_pulse[j] * factor)
                j += 1
                if j == template_pulse.size:
                    break

        pmt_waveform = PMT_Waveform(waveform, pmt_array.get_pmt_object_number(OM_ID))
        pmt_waveform_injected = PMT_Waveform(waveform - injected_data, pmt_array.get_pmt_object_number(OM_ID + 2))
        '''plt.plot(waveform, 'b-')
        plt.plot(waveform - injected_data, 'r-')
        plt.show()'''

        if pmt_waveform.get_pmt_pulse_trigger():
            i_shape = 0
            i_amp = 0
            r_shape = 0
            r_amp = 0

            print(random_time, pmt_waveform_injected.get_pmt_pulse_times())

            for index, value in enumerate(pmt_waveform_injected.get_pmt_pulse_times()):
                if random_time == value:
                    print('found pulse')
                    num[OM_ID] += 1
                    i_shape = pmt_waveform_injected.pmt_waveform_sweep_shape[800 - value]
                    i_amp = pmt_waveform_injected.pmt_waveform_sweep_amp[800 - value]
                    r_shape = pmt_waveform.pmt_waveform_sweep_shape[800 - value]
                    r_amp = pmt_waveform.pmt_waveform_sweep_amp[800 - value]

            # Store random numbers
            t_injected[OM_ID].append(random_time)
            a_injected[OM_ID].append(random_amp)

            injected_array[OM_ID][0].append(i_shape)
            injected_array[OM_ID][1].append(i_amp)
            random_array[OM_ID][0].append(r_shape)
            random_array[OM_ID][1].append(r_amp)

            best_shape = np.amax(pmt_waveform.pmt_waveform_sweep_shape)
            i_best = np.argmax(pmt_waveform.pmt_waveform_sweep_shape)
            best_amplitude = pmt_waveform.pmt_waveform_sweep_amp[i_best]

            x_best_array[OM_ID].append(i_best)
            best_array[OM_ID][0].append(best_shape)
            best_array[OM_ID][1].append(best_amplitude)

            enum[OM_ID] += 1
        iterator += 1
        if iterator == 1000:
            break
        #bar.update(iterator)

    #bar.finish()

    print("Percentage injected found: ",num[0]/enum[0], "% ch0", num[1]/enum[1], "% ch1")

    for i in range(2):
        plt.plot(t_injected[i], injected_array[i][0], 'r.', label="injected")
        plt.plot(t_injected[i], random_array[i][0], 'g.', label="random")
        plt.plot(t_injected[i], injected_array[i][0], 'b.', label="best")
        plt.show()

    output_file.Close()


if __name__ == '__main__':
    main()
