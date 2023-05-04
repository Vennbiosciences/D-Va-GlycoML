import itertools
import os
import random
import re
import time
import warnings
from tkinter import filedialog

from openpyxl import Workbook, load_workbook


### export list of column in txt file
### input (txt file, column header)
### output (column data)
def Func_read_txt_column(txt_filename, column_name):
    txt_file = [line.strip() for line in open(txt_filename).read().splitlines()]
    txt_data = [line.split("\t") for line in txt_file]
    for i in range(len(txt_data[0])):
        if txt_data[0][i] == column_name:
            column_target = i
    return [txt_data[i][column_target] for i in range(1, len(txt_data))]


## program
## program

Input_result = filedialog.askopenfilename(
    initialdir=(os.getcwd()), filetypes=[("pGlyco result", ".txt")], title=("pGlyco result")
)
Input_info = filedialog.askopenfilename(
    initialdir=(os.getcwd()), filetypes=[("MGF info", ".txt")], title=("MGF info")
)

print(Input_result)
print(Input_info)

List_result = Func_read_txt_column(Input_result, "GlySpec")
List_info = Func_read_txt_column(Input_info, "Spectra")

for i in range(len(List_result)):
    List_result[i] = List_result[i].split(".")[0] + "-" + List_result[i].split(".")[1]

for i in range(len(List_info)):
    List_info[i] = List_info[i].split(".")[0] + "-" + List_info[i].split(".")[1]
List_info.insert(0, "Spectra.Spectra")

with open(Input_info.split(".")[0] + "_select.txt", "w") as output:
    line_count = 0
    for line in open(Input_info):
        output_check = False
        if line_count == 0:
            output_check = True
        if List_info[line_count] in List_result:
            output_check = True
        if output_check:
            output.write(line)
        line_count += 1

List_info_selected = Func_read_txt_column(Input_info.split(".")[0] + "_select.txt", "Spectra")
for i in range(len(List_info_selected)):
    List_info_selected[i] = (
        List_info_selected[i].split(".")[0] + "-" + List_info_selected[i].split(".")[1]
    )
List_info_selected.insert(0, "Spectra.Spectra")

with open(Input_info.split(".")[0] + "_select_final.txt", "w") as output:
    line_count = 0
    for line in open(Input_info.split(".")[0] + "_select.txt"):
        output_check = False
        if line_count == 0:
            output_check = True
        if output_check:
            output.write(line)
        line_count += 1
    for i in range(len(List_result)):
        line_count = 0
        for line in open(Input_info.split(".")[0] + "_select.txt"):
            if List_result[i] == List_info_selected[line_count]:
                output.write(line)
            line_count += 1
