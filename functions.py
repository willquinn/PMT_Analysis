import numpy as np


def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser(description="Input file names")
    parser.add_argument('-i', required=True, type=str, help='Input data file')
    parser.add_argument('-config', required=False, type=str, help='Config file')
    args = parser.parse_args()
    return args

def get_normalisation_factor(vector: list):
    norm = 0.0
    for i in range(len(vector)):
        norm += vector[i] * vector[i]
    return np.sqrt(norm)

def get_date_time(input_data_file_name: str):
    temp1 = input_data_file_name.split("/")
    date = temp1[len(temp1) - 2]
    temp2 = temp1[len(temp1) - 1].split("_")
    time = temp2[2]
    return date, time
