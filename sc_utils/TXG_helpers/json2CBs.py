# FOR VISIUM DATA- Extract spot barcodes (CBs) from a Loupe-browser-generated .json file
## Takes in the filepath to the .json (json_path) and the path for the output spot barcode list (CB_path)

## Usage:
## $ path/to/json2CBs.py json_path CB_path whitelist

# def json2CBs(json_path, CB_path, whitelist)
import sys
import json
# from os import path, remove, rename
import pandas as pd

json_path = sys.argv[1]
CB_path = sys.argv[2]
whitelist = sys.argv[3]

# Read in json
with open(json_path) as data_file:
    data = json.load(data_file)

# Read in whitelist, convert long format to wide for easier indexing
coords = pd.read_csv(whitelist, delimiter='\t', names=['CB', 'col_coord', 'row_coord']).pivot(index='row_coord',columns='col_coord', values='CB')

# Get list of row/col coordinates that are included in the tissue section
row_coords = list()
col_coords = list()

for spot in list(data.values())[0]:
    row_coords.append(spot['row'])
    col_coords.append(spot['col'])

# print("Done getting row/col coordinates")

# Pull out CBs from whitelist
CBs_out = list()
for i in list(range(len(row_coords))):
    CBs_out.append(coords.loc()[row_coords[i]+1, col_coords[i]+1]) # whitelist is 1-indexed, json is 0-indexed

# print("Done getting CBs")

# Write CBs as a .txt file
textfile = open(CB_path, "w")
for CB in CBs_out:
    textfile.write(CB + "\n")
textfile.close()
