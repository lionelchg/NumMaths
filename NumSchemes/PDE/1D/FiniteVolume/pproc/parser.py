import numpy as np
import re

def parse_header(header_line):
    p = re.compile('\S+')
    header = []
    while (p.search(header_line)):
        begin, end = p.search(header_line).span()
        header.append(p.search(header_line).group())
        header_line = header_line[end:]
    return header

def read_data(filename, head=False):
    """ Reads data formated in one line header then data.
    The header gives the variables and units in : 
    variable [units] format """
    with open(filename, 'r') as fp:
        header = parse_header(fp.readline().strip())
        data = np.loadtxt(fp)
        if len(data[0, :]) != len(header):
            raise ValueError("Different values of header length and data length")
        indices_sort = np.argsort(data[:, 0])
        for j in range(len(header)):
            data[:, j] = np.array([data[i, j] for i in indices_sort])

    if head:
        return header, data
    else:
        return data