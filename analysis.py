from functions import parse_arguments, process_xml_file, get_date_time
from PMT_Classes import PMT_Array

def main():
    args = parse_arguments()

    input_data_file_name = args.i
    config_file_name = args.config

    date, time = get_date_time(input_data_file_name)

    data_path = "../data/" + date + "/"
    output_path = "../"

    try:
        print(">>> Reading data from file: {}".format(data_path + input_data_file_name))
        date_file = open(data_path + input_data_file_name, 'r')
    except FileNotFoundError as fnf_error:
        print(fnf_error)
        raise Exception("Error opening data file {}".format(data_path + input_data_file_name))

    temp = input_data_file_name.split(".")
    output_file_name = temp[0] + "_output.root"

    try:
        root_file = open(output_path + output_file_name, 'r')
        print(">>> Output file already exists: {}".format(output_path + output_file_name))
    except:
        print(">>> Output doesn't exist, creating {}".format(output_path + output_file_name))

        # Create the object to store all the pmt information
        topology = [2,1]
        pmt_array = PMT_Array(topology, date+"_"+time)
        pmt_array.set_pmt_id("GAO607_"+date+"_"+time, 0)
        pmt_array.set_pmt_id("GAO612_"+date+"_"+time, 1)

        # Set the cuts you wish to apply
        # If you don't do this the defaults are used
        pmt_array.apply_setting(config_file_name)

        process_xml_file(data_path + input_data_file_name, pmt_array)


if __name__ == '__main__':
    main()
