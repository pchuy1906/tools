import numpy as np

import re

import argparse
parser = argparse.ArgumentParser(description='print one line after matching')
# Arguments supported by the code.
parser.add_argument("--file_input",   default='aaa.dat',    help='file input')
parser.add_argument("--keywords",     default='keywords',   help='keywords')

args        = parser.parse_args()
file_input    = args.file_input
keywords      = args.keywords

f = open(file_input, 'rt')

f2 = open("CHNO.dat", "a")


def write_f2(f2, cnames):
    print (cnames)
    longest_cname = max(cnames, key=len)
    print (longest_cname)

    # split the string by "-" and take only the first 4 elements
    list_of_items = longest_cname.split("-")[:4]
    #print (list_of_items)

    numbers = [re.sub(r"\D", "", text) for text in list_of_items]
    #print (numbers)

    string_representation = ' '.join(map(str, numbers))
    print (string_representation)

    f2.write("%s\n" %( string_representation ))


id_collect = 0
while True:
    line = f.readline()
    if line == '':
        write_f2(f2, cnames)
        break

    if keywords in line:
        id_collect += 1
        if id_collect > 1:
            write_f2(f2, cnames)
        cnames = []
    else:
        cname = line.split()[3]
        cnames.append(cname)
