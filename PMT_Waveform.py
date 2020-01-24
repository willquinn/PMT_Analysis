"""
    This class performs all the analysis on a PMT waveform
    Input the relevent PMT and waveform
        - PMT Object ; PMT_Object
        - Waveform ; waveform_array
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from PMT_Object import PMT_Object
import ROOT


class PMT_Waveform:

    def __init__(self, pmt_waveform_list: list, pmt_object: PMT_Object):

        # Store the PMT_Object into memory
        self.pmt_object = pmt_object

        self.pmt_trace_id = self.pmt_object.get_data_id() + "_ch" + self.pmt_object.get_pmt_id() + "_tr" + str(self.pmt_object.get_event_number())

        temp_list_1 = []
        for i in range(len(pmt_waveform_list)):
            if pmt_waveform_list[i] == '' or pmt_waveform_list[i] == '\n':
                pass
            else:
                temp_list_1.append(pmt_waveform_list[i].strip())

        # Store waveform
        self.pmt_waveform = np.array(temp_list_1, dtype='float')

        # Check whether the pulse lies where it is expected i.e above the 100 ns timestamp
        self.pmt_pulse_peak_amplitude = 0.0
        self.pmt_pulse_peak_position = np.argmin(self.pmt_waveform)

        self.pmt_pulse_start = 0
        self.pmt_pulse_end = 0
        self.pmt_baseline = 0.0
        self.pmt_pulse_charge = 0.0
        self.pmt_apulse_charge = 0.0
        self.pmt_pulse_trigger = False
        self.pmt_waveform_reduced = np.array([], dtype='float')

        self.results_dict = {
            "pulse_charge"      : 0.0,
            "pulse_amplitude"   : 0.0,
            "apulse_charge"     : 0.0,
            "mf_pulse_shape"    : 0.0,
            "mf_pulse_amp"      : 0.0
        }

        if self.get_pmt_pulse_peak_position() < self.get_pmt_oject().get_pulse_time_threshold():
            pass
        else:
            # TODO: move this to after cuts. Deal with all cuts here
            self.set_pmt_pulse_trigger(True)

            # Store the baseline
            self.set_pmt_baseline(np.average(self.pmt_waveform[0:self.pmt_pulse_peak_position - 100]))

            # Store the area of the pulse
            # TODO: make this univeral to any size of peak
            self.set_pmt_pulse_start(self.get_pmt_pulse_peak_position() - 20)
            self.set_pmt_pulse_end(self.get_pmt_pulse_peak_position() + 30)

            self.set_pmt_pulse_charge(np.sum(self.get_pmt_waveform()[self.get_pmt_pulse_start():self.get_pmt_pulse_end()] - self.get_pmt_baseline()))
            self.set_pmt_apulse_charge(np.sum(self.get_pmt_waveform()[self.get_pmt_oject().get_apulse_region():] - self.get_pmt_baseline()))
            self.set_pmt_waveform_reduced()
            self.set_pmt_pulse_peak_amplitude(np.amin(self.get_pmt_waveform_reduced()))

            self.results_dict["charge"] = self.get_pmt_pulse_charge()
            self.results_dict["amplitude"] = self.get_pmt_pulse_peak_amplitude()
            self.results_dict["apulse_charge"] = self.get_pmt_apulse_charge()

        # Store the pmt waveform length
        self.pmt_waveform_length = self.pmt_waveform.size
        # self.fill_pmt_hists()

    def get_pulse_trigger(self):
        return self.pmt_pulse_trigger

    def get_pmt_oject(self):
        return self.pmt_object

    def get_pmt_pulse_peak_amplitude(self):
        return self.pmt_pulse_peak_amplitude

    def set_pmt_pulse_peak_amplitude(self, new_amplitude: float):
        self.pmt_pulse_peak_amplitude = new_amplitude

    def get_pmt_apulse_charge(self):
        return self.pmt_apulse_charge

    def set_pmt_apulse_charge(self, new_apulse_charge: float):
        self.pmt_apulse_charge = new_apulse_charge

    def get_pmt_pulse_start(self):
        return self.pmt_pulse_start

    def get_results_dict(self):
        return self.results_dict

    def set_pmt_pulse_start(self, new_pmt_pulse_start: int):
        self.pmt_pulse_start = new_pmt_pulse_start

    def get_pmt_pulse_end(self):
        return self.pmt_pulse_end

    def set_pmt_pulse_end(self, new_pmt_pulse_end: int):
        self.pmt_pulse_end = new_pmt_pulse_end

    def get_pmt_pulse_trigger(self):
        return self.pmt_pulse_trigger

    def set_pmt_pulse_trigger(self, new_pmt_pulse_trigger: bool):
        self.pmt_pulse_trigger = new_pmt_pulse_trigger

    def get_pmt_waveform_reduced(self):
        return self.pmt_waveform_reduced

    def set_pmt_waveform_reduced(self):
        self.pmt_waveform_reduced = self.get_pmt_waveform() - self.get_pmt_baseline()

    def get_pmt_baseline(self):
        return self.pmt_baseline

    def set_pmt_baseline(self, new_pmt_baseline: float):
        self.pmt_baseline = new_pmt_baseline

    def get_pmt_pulse_charge(self):
        return self.pmt_pulse_charge

    def set_pmt_pulse_charge(self, new_pmt_pulse_charge: float):
        self.pmt_pulse_charge = new_pmt_pulse_charge

    def get_pmt_pulse_peak_position(self):
        return self.pmt_pulse_peak_position

    def get_pmt_waveform(self):
        return self.pmt_waveform

    def get_pmt_pulse(self):
        return self.pmt_waveform[self.get_pmt_pulse_start():self.get_pmt_pulse_end()] - self.pmt_baseline

    def get_pmt_waveform_length(self):
        return self.pmt_waveform_length

    def get_pmt_trace_id(self):
        return self.pmt_trace_id

    def fill_pmt_hists(self):
        self.pmt_object.fill_pmt_hists(self.get_results_dict())

    def save_pmt_waveform_histogram(self, root_file: ROOT.TFile):
        name = self.pmt_object.get_pmt_id() + self.get_pmt_trace_id()
        pmt_waveform_hist = ROOT.TH1I(name, name, self.get_pmt_waveform_length(), 0, self.get_pmt_waveform_length())
        for timestamp in range(self.get_pmt_waveform_length()):
            pmt_waveform_hist.SetBinContent(timestamp, self.get_pmt_waveform()[timestamp])
        root_file.cd()
        pmt_waveform_hist.Write()
        del pmt_waveform_hist

    def pmt_pulse_sweep(self, sweep_start: int, sweep_end: int):
        sweep_window_length = self.pmt_object.get_template_pmt_pulse().size

        # plt.plot(self.get_pmt_waveform_reduced())
        # plt.xlabel('timestamp (ns)')
        # plt.title("Waveform {}".format(self.get_pmt_trace_id()))
        # plt.show(block=True)

        matched_filter_shape_list = []
        matched_filter_amplitude_list = []
        for i_sweep in range(sweep_start, sweep_end - sweep_window_length):
            pmt_waveform_section = self.get_pmt_waveform_reduced()[
                                   sweep_start + i_sweep:sweep_start + i_sweep + sweep_window_length]
            if pmt_waveform_section.size == self.pmt_object.get_template_pmt_pulse().size:
                pass
            else:
                print("Section size {} template size {}".format(pmt_waveform_section.size,
                                                                self.pmt_object.get_template_pmt_pulse().size))
                continue

            inner_product = np.dot(self.pmt_object.get_template_pmt_pulse(), pmt_waveform_section)

            matched_filter_amplitude_list.append(inner_product)
            matched_filter_shape_list.append(
                inner_product / np.sqrt(np.dot(pmt_waveform_section, pmt_waveform_section)))

            fig, ax1 = plt.subplots()
            color = 'tab:red'
            ax1.plot(pmt_waveform_section, color=color)
            ax1.set_xlabel('timestamp (ns)')
            ax1.set_ylabel('ADC /mV', color=color)
            ax1.tick_params(axis='y', labelcolor=color)
            color = 'tab:blue'
            ax2 = ax1.twinx()
            ax2.plot(self.pmt_object.get_template_pmt_pulse(), color=color)
            ax2.set_ylabel('Normalised ADC', color=color)  # we already handled the x-label with ax1
            ax2.tick_params(axis='y', labelcolor=color)
            fig.tight_layout()
            plt.text(0, 0,
                     "Shape: {}".format(inner_product / np.sqrt(np.dot(pmt_waveform_section, pmt_waveform_section))))
            plt.show(block=False)
            plt.pause(0.01)
            plt.close()

        matched_filter_shape = np.array(matched_filter_shape_list)
        matched_filter_amplitude = np.array(matched_filter_amplitude_list)

        shape_peaks, _ = find_peaks(matched_filter_shape, height=0.9, distance=int(sweep_window_length / 2))
        amplitude_peaks, _ = find_peaks(matched_filter_amplitude, height=10, distance=int(sweep_window_length / 2))

        fig, ax1 = plt.subplots()

        color = 'tab:red'
        ax1.set_xlabel('timestamp (ns)')
        ax1.set_ylabel('shape', color=color)
        ax1.plot(matched_filter_shape, color=color)
        ax1.plot(shape_peaks, matched_filter_shape[shape_peaks], "x", color='tab:green')
        ax1.tick_params(axis='y', labelcolor=color)
        ax1.plot(np.zeros_like(matched_filter_shape), "--", color="gray")

        ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

        color = 'tab:blue'
        ax2.set_ylabel('Waveform ADC /mV', color=color)  # we already handled the x-label with ax1
        ax2.plot(self.get_pmt_waveform_reduced(), color=color)
        ax2.tick_params(axis='y', labelcolor=color)

        '''ax2.set_ylabel('Amplitude', color=color)  # we already handled the x-label with ax1
        ax2.plot(matched_filter_amplitude, color=color)
        ax2.plot(amplitude_peaks, matched_filter_amplitude[amplitude_peaks], "x", color='tab:orange')
        ax2.tick_params(axis='y', labelcolor=color)'''

        fig.tight_layout()  # otherwise the right y-label is slightly clipped
        plt.show(block=True)


