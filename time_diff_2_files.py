import os
import datetime

fileA="dftb_pin.hsd"
fileB="DFTB_output"
print(datetime.datetime.fromtimestamp(os.path.getmtime(fileB)) - datetime.datetime.fromtimestamp(os.path.getmtime(fileA)))
