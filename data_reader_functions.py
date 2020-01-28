from PMT_Array import PMT_Array
from PMT_Waveform import PMT_Waveform
from xml.dom import minidom
import time as TIME
import ROOT


def process_xml_file(input_data_file_name: str, pmt_array: PMT_Array):
    print(">>> Parsing the data file...")
    processing_start = TIME.time()

    # parse an xml file by name
    xml_file = minidom.parse(input_data_file_name)
    events = xml_file.getElementsByTagName('event')
    parse_time = TIME.time() - processing_start
    print(">>> File is good. Parse time: %.3f s" % parse_time)
    print(">>> Number of Events: {}".format(len(events)))

    save_waveform = True

    event_counter = 0
    percentage_counter = 1
    print(">>> Looping over events")
    temp_start = TIME.time()
    for event_index, event in enumerate(events):

        if event_counter == int(len(events)/10):
            event_counter = 0
            intermediate = TIME.time()
            time_length = intermediate - temp_start
            print(">>>\n>>>  %.3f s.\n" % (intermediate - temp_start))
            temp_start = intermediate
            print("Processed {}% of data...".format(percentage_counter * 10))
            estimate = time_length * (10 - percentage_counter)
            print(">>> Estimated time till termination %.3f seconds\n\n" % estimate)
            percentage_counter += 1

        traces = event.getElementsByTagName('trace')
        for trace_index, trace in enumerate(traces):
            # Channel refers to the pmt number
            channel_id = int(trace.attributes['channel'].value)

            pmt_waveform = PMT_Waveform(trace.firstChild.data.split(" "), pmt_array.get_pmt_object_number(channel_id))

            if save_waveform:
                save_waveform = False
                waveform_root_file = ROOT.TFile("waveform.root", "CREATE")
                pmt_waveform.save_pmt_waveform_histogram(waveform_root_file)
            # check waveform and fill histograms
            if pmt_waveform.get_pulse_trigger():
                pmt_waveform.fill_pmt_hists()

            else:
                pass
            del pmt_waveform

        event_counter += 1


def process_xml_file_new(input_data_file_name: str):
    print(">>> Parsing the data file...")
    processing_start = TIME.time()

    # parse an xml file by name
    xml_file = minidom.parse(input_data_file_name)
    events = xml_file.getElementsByTagName('event')
    parse_time = TIME.time() - processing_start
    intermediate_time = TIME
    print(">>> File is good. Parse time: %.3f s" % parse_time)
    print(">>> Number of Events: {}".format(len(events)))

    waveform_data_list = []
    for event_index, event in enumerate(events):

        traces = event.getElementsByTagName('trace')
        for trace_index, trace in enumerate(traces):
            waveform_data_list.append([int(trace.attributes['channel'].value), trace.firstChild.data.split(" ")])
    return waveform_data_list


def read_test_file():
    try:
        pmt_data_file = open("./test_data.dat", 'r')
    except FileNotFoundError as fnf_error:
        print(fnf_error)
        raise Exception("Error opening data file.")
    data_list = []
    for line_index, line in enumerate(pmt_data_file.readlines()):
        waveform = []
        for token_index, token in enumerate(line.split()):
            waveform.append(token.strip())
        data_list.append(waveform)
    return data_list





