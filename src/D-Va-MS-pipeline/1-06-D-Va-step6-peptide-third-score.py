import itertools
import math
import os
import random
import time
import warnings
from tkinter import filedialog

import matplotlib.pyplot as Spe_plt
from openpyxl import Workbook, load_workbook

warnings.filterwarnings("ignore")

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


### median
### input (list of numbers)
### output (median value)
def Func_median(numbers):
    numbers = sorted(numbers)
    center = len(numbers) / 2
    if len(numbers) % 2 == 0:
        return sum(numbers[center - 1 : center + 1]) / 2.0
    else:
        return numbers[center]


### third
### input (list of numbers)
### output (third value)
def Func_third(numbers):
    numbers = sorted(numbers)
    third = len(numbers) // 3
    return numbers[third]


### nonredundant list
### input (data list)
### output (nonredundant list)
def Func_delete_redundant(input):
    input_data = input
    output_data = []
    for i in range(len(input_data)):
        if input_data[i] not in output_data:
            output_data.append(input_data[i])
    return output_data


### Column head of glycan
### input (txt file)
### output (glycan head)
def Func_extract_glycan_head(file):
    Extract_file = [line.strip() for line in open(file).read().splitlines()]
    for i in range(len(Extract_file)):
        Extract_file[i] = Extract_file[i].split("\t")
    for i in range(len(Extract_file[0])):
        if "Glycan(" in Extract_file[0][i] and "Corrected" not in Extract_file[0][i]:
            glycan_head = Extract_file[0][i]
    return glycan_head


## program
## program

# input-directory, txt files
Input_filename = filedialog.askopenfilename(
    initialdir=(os.getcwd()), filetypes=[("input data", ".txt")], title=("input data")
)
print(Input_filename)

Input_peptide = Func_read_txt_column(Input_filename, "Peptide")
Input_glycan = Func_read_txt_column(Input_filename, Func_extract_glycan_head(Input_filename))
Input_score = Func_read_txt_column(Input_filename, "TotalScore")

List_peptide = Func_delete_redundant(Input_peptide)
List_peptide_median_score = []

# output third-highest glycopeptide score for each peptide backbone
for i in range(len(List_peptide)):
    Current_score_pool = []
    for j in range(len(Input_peptide)):
        if List_peptide[i] == Input_peptide[j]:
            Current_score_pool.append(float(Input_score[j]))
    List_peptide_median_score.append(Func_third(Current_score_pool))

List_GP_final = []

# output each glycopeptide pass step6
with open(Input_filename.split(".")[0] + "-step6-output-01-check.txt", "w") as Output:
    Output.write(
        "Peptide\t" + Func_extract_glycan_head(Input_filename) + "\tTotalScore\tD-Va-Step6-check\n"
    )
    for i in range(len(Input_peptide)):
        Output_value = 1
        for j in range(len(List_peptide)):
            if (
                Input_peptide[i] == List_peptide[j]
                and float(Input_score[i]) < List_peptide_median_score[j]
            ):
                Output_value = 0
        if Output_value == 1:
            List_GP_final.append(str(Input_peptide[i]) + " - " + str(Input_glycan[i]))
        Output.write(
            str(Input_peptide[i])
            + "\t"
            + str(Input_glycan[i])
            + "\t"
            + str(Input_score[i])
            + "\t"
            + str(Output_value)
        )
        Output.write("\n")

# output unique glycopeptide list
List_GP_final = Func_delete_redundant(List_GP_final)
with open(Input_filename.split(".")[0] + "-step6-output-02-GP-list.txt", "w") as Output:
    Output.write("Glycopeptide\n")
    for i in range(len(List_GP_final)):
        Output.write(List_GP_final[i])
        Output.write("\n")

# output final PSM passed step6
with open(Input_filename.split(".")[0] + "-step6-output-03-psm.txt", "w") as Output:
    Extract_file = [line.strip() for line in open(Input_filename).read().splitlines()]
    Output.write(Extract_file[0])
    Output.write("\n")
    for i in range(1, len(Extract_file)):
        Output_check = False
        for j in range(len(List_GP_final)):
            if (
                Input_peptide[i - 1] == List_GP_final[j].split(" - ")[0]
                and Input_glycan[i - 1] == List_GP_final[j].split(" - ")[1]
            ):
                Output_check = True
        if Output_check:
            Output.write(Extract_file[i])
            Output.write("\n")
