from PMT_Classes import PMT_Object, PMT_Waveform, PMT_Array
import ROOT
import numpy as np


def main():
    #create_template_histogram_temp()
    run_number = "214"
    data_file_name = "/Users/willquinn/Documents/PhD/SuperNEMO/SNEMO_ComData_Analysis/GV_XW_ComData/Data/run_" + run_number + "/run_" + run_number
    template_file_name = "/Users/willquinn/Documents/PhD/SuperNEMO/SNEMO_ComData_Analysis/template.root"
    topology = [20, 20]

    pmt_array = PMT_Array(topology, run_number)
    pmt_array.set_pmt_templates(template_file_name)

    Read_Data(data_file_name, topology, pmt_array)


def Read_Data(pmt_data_filename: str, topology: list, pmt_array: PMT_Array):

    try:
        pmt_data_file = open(pmt_data_filename, 'r')
    except FileNotFoundError as fnf_error:
        print(fnf_error)
        raise Exception("Error opening data file. Skip to the next file...")

    new_waveform_bool = False
    line_number_int = 0
    waveform_number_int = 0

    for pmt_data_index, pmt_data_line in enumerate(pmt_data_file.readlines()[10:]):
        pmt_data_line_tokens = pmt_data_line.split(" ")

        if pmt_data_line_tokens[0] == "=" and pmt_data_line_tokens[1] == "HIT":
            new_waveform_bool = True
            line_number_int = 0
        else:
            pass

        if new_waveform_bool and line_number_int == 1:
            pmt_slot_number = int(pmt_data_line_tokens[1])
            pmt_channel_number = int(pmt_data_line_tokens[3])
            pmt_number = int(pmt_slot_number*topology[0] + pmt_channel_number)
            event_id_LTO = int(pmt_data_line_tokens[5])
            event_id_HT = int(pmt_data_line_tokens[7])
            pmt_waveform_peak_cell = int(pmt_data_line_tokens[27])
            pmt_waveform_charge = float(pmt_data_line_tokens[29])
            pmt_waveform_rise_time = float(pmt_data_line_tokens[39])
            if int(event_id_HT) != 0:
                pass
            elif int(event_id_LTO) != 0:
                pass
            else:
                pass

        elif new_waveform_bool and line_number_int == 2:
            if event_id_HT != 0:
                pmt_adc_values = []
                for i_adc in range(len(pmt_data_line_tokens)):
                    pmt_adc_values.append(pmt_data_line_tokens[i_adc])

                pmt_waveform = PMT_Waveform(pmt_adc_values, pmt_array.get_pmt_object_position([pmt_slot_number, pmt_channel_number]))
                pmt_waveform.fill_pmt_hists()
                pmt_waveform.pmt_pulse_sweep(0, 500)

                del pmt_waveform

            new_waveform_bool = False

        line_number_int += 1



    pmt_data_file.close()


def create_template_histogram_temp():
    source_file = ROOT.TFile("/Users/willquinn/Documents/PhD/SuperNEMO/SNEMO_ComData_Analysis/GV_XW_ComData/Output_files/Average_Waveforms_run214.root", "READ")
    source_waveform = ROOT.TH1F(source_file.Get("Waveform_10_3_0"))
    #source_file.Close()

    output_file = ROOT.TFile("/Users/willquinn/Documents/PhD/SuperNEMO/SNEMO_ComData_Analysis/template.root", "RECREATE")
    template_hist = ROOT.TH1F("template", "template", 80, 0, 80)
    output_file.cd()

    pmt_object = PMT_Object("template", "template")

    temp_list = []
    for i_bin in range(int(source_waveform.GetEntries())):
        temp_list.append(str(source_waveform.GetBinContent(i_bin)))
        #print(i_bin, " ", str(source_waveform.GetBinContent(i_bin + 1)))
    pmt_waveform = PMT_Waveform(temp_list, pmt_object)

    for i_bin in range(80):
        template_hist.SetBinContent(i_bin+1, pmt_waveform.get_pmt_waveform()[260 + i_bin] - pmt_waveform.get_pmt_baseline())
    output_file.Write()
    output_file.Close()

if __name__ == '__main__':
    main()
