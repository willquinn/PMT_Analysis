import ROOT
from PMT_Classes import PMT_Waveform, PMT_Object, PMT_Array
from xml.dom import minidom
import time as TIME
import numpy as np


def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser(description="Input file names")
    parser.add_argument('-i', required=True, type=str, help='Input data file')
    parser.add_argument('-config', required=False, type=str, help='Config file')
    args = parser.parse_args()
    return args


def get_date_time(input_data_file_name: str):
    temp1 = input_data_file_name.split("/")
    date = temp1[len(temp1) - 2]
    temp2 = temp1[len(temp1) - 1].split("_")
    time = temp2[2]
    return date, time


def process_xml_file(input_data_file_name: str, pmt_array: PMT_Array):
    print(">>> Parsing the data file...")
    processing_start = TIME.time()

    # parse an xml file by name
    xml_file = minidom.parse(input_data_file_name)
    events = xml_file.getElementsByTagName('event')
    parse_time = TIME.time() - processing_start
    print(">>> File is good. Parse time: %.3f s" % parse_time)
    print(">>> Number of Events: {}".format(len(events)))

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

            # check waveform and fill histograms
            if pmt_waveform.get_pulse_trigger():
                pmt_waveform.fill_pmt_hists()
            else:
                pass
            del pmt_waveform

        event_counter += 1

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