import sys
sys.path.insert(1, '..')

from scr.PMT_Array import PMT_Array
from functions.other_functions import sncd_parse_arguments, get_run_number, get_data_path
from functions.data_reader_functions import process_crd_file, read_filenames
import time as TIME
import ROOT


def main():
    args = sncd_parse_arguments()
    data_file_names = args.i
    config_file_name = args.c
    sweep_bool = args.sweep
    template_file_names = args.t
    output_file_name = args.o

    run_number = get_run_number(data_file_names)

    topology = [20, 13]

    template_list = read_filenames(template_file_names)
    list_of_data_file_names = read_filenames(data_file_names)

    pmt_array = PMT_Array(topology, run_number)
    if config_file_name is not None:
        pmt_array.apply_setting(config_file_name)
    if template_file_names is not None:
        pmt_array.set_pmt_templates(template_list, "template")
    if sweep_bool == "True":
        pmt_array.set_sweep_bool(True)

    print(">>> Looping of files...")

    data_path = get_data_path(data_file_names)

    waveform_output_root_file = ROOT.TFile("waveforms.root", "RECREATE")

    temp_start = TIME.time()
    for index, data_file_name in enumerate(list_of_data_file_names):
        process_crd_file(data_path + data_file_name, pmt_array, waveform_output_root_file)

        intermediate = TIME.time()
        time_length = intermediate - temp_start
        print(">>>\n>>>  %.3f s.\n" % (intermediate - temp_start))
        temp_start = intermediate
        print("Processed {} of the {} data files...".format(index + 1, len(list_of_data_file_names)))
        estimate = time_length * (len(list_of_data_file_names) - index - 1)
        print(">>> Estimated time till termination %.3f seconds\n\n" % estimate)

    pmt_array.save_to_file(output_file_name)


def create_template_histogram_temp():
    source_file = ROOT.TFile(
        "/Users/willquinn/Documents/PhD/SuperNEMO/SNEMO_ComData_Analysis/GV_XW_ComData/Output_files/Average_Waveforms_run214.root",
        "READ")
    source_waveform = ROOT.TH1F(source_file.Get("Waveform_10_3_0"))
    # source_file.Close()

    output_file = ROOT.TFile("/Users/willquinn/Documents/PhD/SuperNEMO/SNEMO_ComData_Analysis/template.root",
                             "RECREATE")
    template_hist = ROOT.TH1F("template", "template", 80, 0, 80)
    output_file.cd()

    pmt_object = PMT_Object("template", "template")

    temp_list = []
    for i_bin in range(int(source_waveform.GetEntries())):
        temp_list.append(str(source_waveform.GetBinContent(i_bin)))
        # print(i_bin, " ", str(source_waveform.GetBinContent(i_bin + 1)))
    pmt_waveform = PMT_Waveform(temp_list, pmt_object)

    for i_bin in range(80):
        template_hist.SetBinContent(i_bin + 1,
                                    pmt_waveform.get_pmt_waveform()[260 + i_bin] - pmt_waveform.get_pmt_baseline())
    output_file.Write()
    output_file.Close()


if __name__ == '__main__':
    main()
