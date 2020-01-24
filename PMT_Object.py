"""
    This is the base class for a Photomultiplier Tube
    To create an unique object:
        - PMT identifier ; PMT_id
        - Data identifier ; Data_id
    This will contain the results of waveform analysis from class PMT_Waveform
"""

import ROOT


class PMT_Object:

    def __init__(self, pmt_id: str, data_id: str):
        self.pmt_id = pmt_id
        self.data_id = data_id

        self.number_of_events = 0

        self.charge_cut = 6  # pC
        self.sweep_range = [0, 500]
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
        self.charge_range = [0, 100]
        self.amp_range = [0, 100]
        self.mf_shape_range = [-1, 1]
        self.mf_amp_range = [0, 100]
        self.apulse_region = 500

        self.pmt_pulse_charge_hist = ROOT.TH1F(self.get_pmt_id() + "_pulse_charge_spectrum",
                                               self.get_pmt_id() + "_pulse_charge_spectrum", self.nbins,
                                               self.charge_range[0], self.charge_range[1])
        self.pmt_pulse_amplitude_hist = ROOT.TH1F(self.get_pmt_id() + "_pulse_amplitude_spectrum",
                                                  self.get_pmt_id() + "_pulse_amplitude_spectrum", self.nbins,
                                                  self.amp_range[0], self.amp_range[1])
        self.pmt_apulse_charge_hist = ROOT.TH1F(self.get_pmt_id() + "_apulse_charge_spectrum",
                                                self.get_pmt_id() + "_apulse_charge_spectrum", self.nbins,
                                                self.charge_range[0], self.charge_range[1])
        self.pmt_pulse_mf_shape_hist = ROOT.TH1F(self.get_pmt_id() + "_mf_shape", self.get_pmt_id() + "_mf_shape",
                                                 self.nbins, self.mf_shape_range[0], self.mf_shape_range[1])
        self.pmt_pulse_mf_amplitude_hist = ROOT.TH1F(self.get_pmt_id() + "_mf_amplitude",
                                                     self.get_pmt_id() + "_mf_amplitude", self.nbins,
                                                     self.mf_amp_range[0], self.mf_amp_range[1])
        self.pmt_pulse_mf_shape_mf_amplitude_hist = ROOT.TH2F(self.get_pmt_id() + "_mf_shape_mf_amplitude",
                                                              self.get_pmt_id() + "_mf_shape_mf_amplitude", self.nbins,
                                                              self.mf_shape_range[0], self.mf_shape_range[1],
                                                              self.nbins, self.mf_amp_range[0], self.mf_amp_range[1])
        self.pmt_pulse_mf_shape_p_amplitude_hist = ROOT.TH2F(self.get_pmt_id() + "_mf_shape_p_amplitude",
                                                             self.get_pmt_id() + "_mf_shape_p_amplitude", self.nbins,
                                                             self.mf_shape_range[0], self.mf_shape_range[1], self.nbins,
                                                             self.amp_range[0], self.amp_range[1])
        self.pmt_pulse_mf_amplitude_p_amplitude_hist = ROOT.TH2F(self.get_pmt_id() + "_mf_amplitude_p_amplitude",
                                                                 self.get_pmt_id() + "_mf_amplitude_p_amplitude",
                                                                 self.nbins, self.mf_amp_range[0], self.mf_amp_range[1],
                                                                 self.nbins, self.amp_range[0], self.amp_range[1])

        self.histogram_dict = {
            "pulse_charge_spectrum": self.pmt_pulse_charge_hist,
            "pulse_amplitude_spectrum": self.pmt_pulse_amplitude_hist,
            "apulse_charge_spectrum": self.pmt_apulse_charge_hist,
            "mf_shape": self.pmt_pulse_mf_shape_hist,
            "mf_amplitude": self.pmt_pulse_mf_amplitude_hist,
            "mf_shape_mf_amplitude": self.pmt_pulse_mf_shape_mf_amplitude_hist,
            "mf_shape_p_amplitude": self.pmt_pulse_mf_shape_p_amplitude_hist,
            "mf_amplitude_p_amplitude": self.pmt_pulse_mf_amplitude_p_amplitude_hist
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
                                               self.get_pmt_id() + "_pulse_charge_spectrum", self.nbins,
                                               self.charge_range[0],
                                               self.charge_range[1])
        self.pmt_pulse_amplitude_hist = ROOT.TH1F(self.get_pmt_id() + "_pulse_amplitude_spectrum",
                                                  self.get_pmt_id() + "_pulse_amplitude_spectrum", self.nbins,
                                                  self.amp_range[0], self.amp_range[1])
        self.pmt_apulse_charge_hist = ROOT.TH1F(self.get_pmt_id() + "_apulse_charge_spectrum",
                                                self.get_pmt_id() + "_apulse_charge_spectrum", self.nbins,
                                                self.charge_range[0], self.charge_range[1])
        self.pmt_pulse_mf_shape_hist = ROOT.TH1F(self.get_pmt_id() + "_mf_shape", self.get_pmt_id() + "_mf_shape",
                                                 self.nbins,
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
                                                                 self.get_pmt_id() + "_mf_amplitude_p_amplitude",
                                                                 self.nbins,
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
