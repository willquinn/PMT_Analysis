import ROOT
from PMT_Object import PMT_Object


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
