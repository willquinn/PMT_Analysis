import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import ROOT


def get_normalisation_factor(vector: list):
    norm = 0.0
    for i in range(len(vector)):
        norm += vector[i] * vector[i]
    return np.sqrt(norm)


class PMT_Object:
    """
    This class is for a Photomultiplier Tube
    To create an unique object:
        - PMT identifier ; PMT_id
        - Data identifier ; Data_id
    This will contain the results of waveform analysis form class PMT_Waveform
    """

    def __init__(self, pmt_id: str, data_id: str):
        self.pmt_id = pmt_id
        self.data_id = data_id

        self.number_of_events = 0

        self.charge_cut = 6  # pC
        self.sweep_range = [0,500]
        self.pulse_time_threshold = 100
        self.template_pmt_pulse = np.array([], dtype='float')

        self.histogram_names_list = ["pulse_charge_spectrum",
                                     "pulse_amplitude_spectrum",
                                     "apulse_charge_spectrum",
                                     "mf_shape",
                                     "mf_amplitude",
                                     "mf_shape_mf_amplitude",
                                     "mf_shape_p_amplitude",
                                     "mf_amplitude_p_amplitude"]
        
        self.nbins = 100
        self.charge_range = [0,100]
        self.amp_range = [0,100]
        self.mf_shape_range = [-1,1]
        self.mf_amp_range = [0,100]
        self.apulse_region = 500

        self.pmt_pulse_charge_hist = ROOT.TH1F(self.get_pmt_id() + "_pulse_charge_spectrum", self.get_pmt_id() + "_pulse_charge_spectrum", self.nbins, self.charge_range[0], self.charge_range[1])
        self.pmt_pulse_amplitude_hist = ROOT.TH1F(self.get_pmt_id() + "_pulse_amplitude_spectrum", self.get_pmt_id() + "_pulse_amplitude_spectrum", self.nbins, self.amp_range[0], self.amp_range[1])
        self.pmt_apulse_charge_hist = ROOT.TH1F(self.get_pmt_id() + "_apulse_charge_spectrum", self.get_pmt_id() + "_apulse_charge_spectrum",self.nbins, self.charge_range[0], self.charge_range[1])
        self.pmt_pulse_mf_shape_hist = ROOT.TH1F(self.get_pmt_id() + "_mf_shape", self.get_pmt_id() + "_mf_shape", self.nbins, self.mf_shape_range[0], self.mf_shape_range[1])
        self.pmt_pulse_mf_amplitude_hist = ROOT.TH1F(self.get_pmt_id() + "_mf_amplitude", self.get_pmt_id() + "_mf_amplitude", self.nbins, self.mf_amp_range[0], self.mf_amp_range[1])
        self.pmt_pulse_mf_shape_mf_amplitude_hist = ROOT.TH2F(self.get_pmt_id() + "_mf_shape_mf_amplitude", self.get_pmt_id() + "_mf_shape_mf_amplitude", self.nbins, self.mf_shape_range[0], self.mf_shape_range[1], self.nbins, self.mf_amp_range[0], self.mf_amp_range[1])
        self.pmt_pulse_mf_shape_p_amplitude_hist = ROOT.TH2F(self.get_pmt_id() + "_mf_shape_p_amplitude", self.get_pmt_id() + "_mf_shape_p_amplitude", self.nbins, self.mf_shape_range[0], self.mf_shape_range[1], self.nbins, self.amp_range[0], self.amp_range[1])
        self.pmt_pulse_mf_amplitude_p_amplitude_hist = ROOT.TH2F(self.get_pmt_id() + "_mf_amplitude_p_amplitude", self.get_pmt_id() + "_mf_amplitude_p_amplitude", self.nbins, self.mf_amp_range[0], self.mf_amp_range[1], self.nbins, self.amp_range[0], self.amp_range[1])

        self.histogram_dict = {
            "pulse_charge_spectrum"     : self.pmt_pulse_charge_hist,
            "pulse_amplitude_spectrum"  : self.pmt_pulse_amplitude_hist,
            "apulse_charge_spectrum"    : self.pmt_apulse_charge_hist,
            "mf_shape"                  : self.pmt_pulse_mf_shape_hist,
            "mf_amplitude"              : self.pmt_pulse_mf_amplitude_hist,
            "mf_shape_mf_amplitude"     : self.pmt_pulse_mf_shape_mf_amplitude_hist,
            "mf_shape_p_amplitude"      : self.pmt_pulse_mf_shape_p_amplitude_hist,
            "mf_amplitude_p_amplitude"  : self.pmt_pulse_mf_amplitude_p_amplitude_hist
        }

    def get_apulse_region(self):
        return self.apulse_region

    def set_apulse_region(self, new_apulse_pos: int):
        self.apulse_region = new_apulse_pos

    def get_histogram_dict(self):
        return self.histogram_dict

    def get_histogram(self, description: str):
        return self.get_histogram_dict()[description]

    def get_sweep_range(self):
        return self.sweep_range

    def get_sweep_range_min(self):
        return self.sweep_range[0]

    def get_sweep_range_max(self):
        return self.sweep_range[1]

    def set_sweep_range(self, new_sweep_range: list):
        assert len(new_sweep_range) == 2
        assert new_sweep_range[1] > new_sweep_range[0]

    def get_pulse_time_threshold(self):
        return self.pulse_time_threshold

    def set_pulse_time_threshold(self, new_threshold: int):
        self.pulse_time_threshold = new_threshold

    def update_hists(self):
        self.pmt_pulse_charge_hist = ROOT.TH1F(self.get_pmt_id() + "_pulse_charge_spectrum",
                                               self.get_pmt_id() + "_pulse_charge_spectrum", self.nbins, self.charge_range[0],
                                               self.charge_range[1])
        self.pmt_pulse_amplitude_hist = ROOT.TH1F(self.get_pmt_id() + "_pulse_amplitude_spectrum",
                                                  self.get_pmt_id() + "_pulse_amplitude_spectrum", self.nbins,
                                                  self.amp_range[0], self.amp_range[1])
        self.pmt_apulse_charge_hist = ROOT.TH1F(self.get_pmt_id() + "_apulse_charge_spectrum",
                                                self.get_pmt_id() + "_apulse_charge_spectrum", self.nbins,
                                                self.charge_range[0], self.charge_range[1])
        self.pmt_pulse_mf_shape_hist = ROOT.TH1F(self.get_pmt_id() + "_mf_shape", self.get_pmt_id() + "_mf_shape", self.nbins,
                                                 self.mf_shape_range[0], self.mf_shape_range[1])
        self.pmt_pulse_mf_amplitude_hist = ROOT.TH1F(self.get_pmt_id() + "_mf_amplitude", self.pmt_id + "_mf_amplitude",
                                                     self.nbins, self.mf_amp_range[0], self.mf_amp_range[1])
        self.pmt_pulse_mf_shape_mf_amplitude_hist = ROOT.TH2F(self.get_pmt_id() + "_mf_shape_mf_amplitude",
                                                              self.get_pmt_id() + "_mf_shape_mf_amplitude", self.nbins,
                                                              self.mf_shape_range[0], self.mf_shape_range[1],
                                                              self.nbins, self.mf_amp_range[0], self.mf_amp_range[1])
        self.pmt_pulse_mf_shape_p_amplitude_hist = ROOT.TH2F(self.get_pmt_id() + "_mf_shape_p_amplitude",
                                                             self.get_pmt_id() + "_mf_shape_p_amplitude", self.nbins,
                                                             self.mf_shape_range[0], self.mf_shape_range[1], self.nbins,
                                                             self.amp_range[0], self.amp_range[1])
        self.pmt_pulse_mf_amplitude_p_amplitude_hist = ROOT.TH2F(self.get_pmt_id() + "_mf_amplitude_p_amplitude",
                                                                 self.get_pmt_id() + "_mf_amplitude_p_amplitude", self.nbins,
                                                                 self.mf_amp_range[0], self.mf_amp_range[1], self.nbins,
                                                                 self.amp_range[0], self.amp_range[1])

    def get_nbins(self):
        return self.nbins

    def set_nbins(self, new_nbins):
        self.nbins = new_nbins

    def get_charge_range(self):
        return self.charge_range

    def set_charge_range(self, new_charge_range: list):
        assert len(new_charge_range) == 2
        assert new_charge_range[1] > new_charge_range[0]
        self.charge_range = new_charge_range

    def get_amp_range(self):
        return self.amp_range

    def set_amp_range(self, new_amp_range: list):
        assert len(new_amp_range) == 2
        assert new_amp_range[1] > new_amp_range[0]
        self.amp_range = new_amp_range

    def get_mf_shape_range(self):
        return self.mf_shape_range

    def set_mf_shape_range(self, new_mf_shape_range: list):
        assert len(new_mf_shape_range) == 2
        assert new_mf_shape_range[1] > new_mf_shape_range[0]
        self.mf_shape_range = new_mf_shape_range

    def get_mf_amp_range(self):
        return self.mf_amp_range

    def set_mf_amp_range(self, new_mf_amp_range: list):
        assert len(new_mf_amp_range) == 2
        assert new_mf_amp_range[1] > new_mf_amp_range[0]
        self.mf_amp_range = new_mf_amp_range

    def get_histogram_names(self):
        return self.histogram_names_list

    def get_histogram_names_dict(self):
        return self.histogram_dict

    def get_pmt_id(self):
        return self.pmt_id

    def set_pmt_id(self, new_pmt_id: str):
        self.pmt_id = new_pmt_id

    def get_data_id(self):
        return self.data_id

    def get_charge_cut(self):
        return self.charge_cut

    def get_event_number(self):
        return self.number_of_events

    def set_number_of_events(self, new_event_number):
        self.number_of_events = new_event_number

    def create_pmt_pulse_template(self, root_file_name: str, template_histogram_name: str):
        template_root_file = ROOT.TFile(root_file_name, "READ")
        template_histogram = ROOT.TH1F(template_root_file.Get(template_histogram_name))
        template_list = []
        for i_bin in range(int(template_histogram.GetEntries())):
            template_list.append(template_histogram.GetBinContent(i_bin))

        norm = get_normalisation_factor(template_list)

        self.set_template_pmt_pulse(np.array(template_list, dtype='float') / norm)
        # print(self.get_template_pmt_pulse())

    def set_template_pmt_pulse(self, new_template_pmt_pulse: np.array):
        self.template_pmt_pulse = new_template_pmt_pulse

    def get_template_pmt_pulse(self):
        return self.template_pmt_pulse

    def get_pmt_pulse_charge_hist(self):
        return self.pmt_pulse_charge_hist

    def fill_pmt_pulse_charge_hist(self, value: float):
        self.pmt_pulse_charge_hist.Fill(value)

    def get_pmt_pulse_amplitude_hist(self):
        return self.pmt_pulse_amplitude_hist

    def fill_pmt_pulse_amplitude_hist(self, value: float):
        self.pmt_pulse_amplitude_hist.Fill(value)

    def get_pmt_pulse_mf_shape_hist(self):
        return self.pmt_pulse_mf_shape_hist

    def fill_pmt_pulse_mf_shape_hist(self, value: float):
        self.pmt_pulse_mf_shape_hist.Fill(value)

    def get_pmt_pulse_mf_amplitude_hist(self):
        return self.pmt_pulse_mf_amplitude_hist

    def fill_pmt_pulse_mf_amplitude_hist(self, value: float):
        self.pmt_pulse_mf_amplitude_hist.Fill(value)

    def get_pmt_pulse_mf_shape_mf_amplitude(self):
        return self.pmt_pulse_mf_shape_mf_amplitude_hist

    def fill_pmt_pulse_mf_shape_mf_amplitude(self, x_value: float, y_value: float):
        self.pmt_pulse_mf_shape_mf_amplitude_hist.Fill(x_value, y_value)

    def get_pmt_pulse_mf_shape_p_amplitude_hist(self):
        return self.pmt_pulse_mf_shape_p_amplitude_hist

    def fill_pmt_pulse_mf_shape_p_amplitude_hist(self, x_value: float, y_value: float):
        self.pmt_pulse_mf_shape_p_amplitude_hist.Fill(x_value, y_value)

    def get_pmt_pulse_mf_amplitude_p_amplitude_hist(self):
        return self.pmt_pulse_mf_amplitude_p_amplitude_hist

    def fill_pmt_pulse_mf_amplitude_p_amplitude_hist(self, x_value, y_value):
        self.pmt_pulse_mf_amplitude_p_amplitude_hist.Fill(x_value, y_value)

    def get_pmt_apulse_charge_hist(self):
        return self.pmt_apulse_charge_hist

    def fill_pmt_apulse_charge_hist(self, value: float):
        self.pmt_apulse_charge_hist.Fill(value)

    def fill_pmt_hists(self, results: dict):

        pulse_charge: float = results["pulse_charge"]
        pulse_amplitude: float = results["pulse_amplitude"]
        apulse_charge: float = results["apulse_charge"]
        mf_amplitude: float = results["mf_pulse_amp"]
        mf_shape: float = results["mf_pulse_shape"]

        self.fill_pmt_pulse_charge_hist(pulse_charge)
        self.fill_pmt_pulse_amplitude_hist(pulse_amplitude)
        self.fill_pmt_apulse_charge_hist(apulse_charge)
        self.fill_pmt_pulse_mf_shape_hist(mf_shape)
        self.fill_pmt_pulse_mf_amplitude_hist(mf_amplitude)
        self.fill_pmt_pulse_mf_shape_mf_amplitude(mf_shape, mf_amplitude)
        self.fill_pmt_pulse_mf_shape_p_amplitude_hist(pulse_amplitude, mf_shape)
        self.fill_pmt_pulse_mf_amplitude_p_amplitude_hist(pulse_amplitude, mf_amplitude)

    def save_to_file(self, root_file_name: str):
        root_file = ROOT.TFile(root_file_name, "RECREATE")
        root_file.cd()
        for hist in self.get_histogram_dict().keys():
            self.get_histogram_dict()[hist].Write()
        root_file.Close()

    def save_histograms(self, directory: ROOT.TDirectory):
        directory.cd()
        for hist in self.get_histogram_dict().keys():
            self.get_histogram_dict()[hist].Write()

    def save_histogram(self, root_file: ROOT.TFile, hist, write_function: str):
        if write_function in ['RECREATE', "CREATE", "UPDATE"]:
            pass
        else:
            print("Invalid write function.")

class PMT_Waveform:
    """
    This class performs all the analysis on a PMT waveform
    Input the relevent PMT and waveform
        - PMT Object ; PMT_Object
        - Waveform ; waveform_array
    """

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


class PMT_Array:

    def __init__(self, topology: list, data_id: str):
        self.pmt_topology = topology
        self.pmt_object_array = []
        self.pmt_total_number = topology[0] * topology[1]

        for i in range(topology[0]):
            for j in range(topology[1]):
                pmt_number = i * topology[0] + j
                pmt_object = PMT_Object(str(pmt_number), data_id)
                self.append_pmt_object_array(pmt_object)
                del pmt_object

    def save_to_file(self, output_root_filename: str):
        output_root_file = ROOT.TFile(output_root_filename, "RECREATE")

        for i_pmt in range(self.get_pmt_total_number()):
            output_root_file.cd()
            output_root_file.mkdir(self.get_pmt_object_number(i_pmt).get_pmt_id())
            directory = output_root_file.GetDirectory(self.get_pmt_object_number(i_pmt).get_pmt_id())
            self.get_pmt_object_number(i_pmt).save_histograms(directory)

    def get_pmt_topology(self):
        return self.pmt_topology

    def get_pmt_oject_array(self):
        return self.pmt_object_array

    def get_pmt_total_number(self):
        return self.pmt_total_number

    def append_pmt_object_array(self, pmt_object: PMT_Object):
        self.pmt_object_array.append(pmt_object)

    def get_pmt_object_position(self, pmt_position: list):
        assert len(pmt_position) < 3
        if len(pmt_position) == 1:
            pmt_number = pmt_position[0]
        else:
            pmt_number = pmt_position[0] * self.get_pmt_topology()[0] + pmt_position[1]
        return self.pmt_object_array[pmt_number]

    def get_pmt_object_number(self, pmt_number: int):
        return self.get_pmt_oject_array()[pmt_number]

    def set_pmt_templates(self, template_root_file_name: str):
        for i_pmt in range(self.get_pmt_total_number()):
            self.get_pmt_object_number(i_pmt).create_pmt_pulse_template(template_root_file_name, template_root_file_name)

    def apply_setting(self, config_file_name: str):
        try:
            config_file = open(config_file_name,'r')
        except FileNotFoundError as fnf_error:
            print(fnf_error)
            raise Exception("Error opening config file")


    def set_pmt_id(self, pmt_id: str, pmt_object_number: int):
        self.get_pmt_object_number(pmt_object_number).set_pmt_id(pmt_id)
