import math
import os
import sys
from tkinter import filedialog


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


## program
## program

# input-txt files of mzML scan info
Input_filename = filedialog.askopenfilename(
    initialdir=(os.getcwd()), filetypes=[("input data", ".txt")], title=("input data")
)
print(Input_filename)

List_file = Func_read_txt_column(Input_filename, "File")
List_file_nonredundant = Func_delete_redundant(List_file)

# output_mzML scan info count
with open(Input_filename.split(".")[0] + "-count.txt", "w") as output_result:
    for i in range(len(List_file_nonredundant)):
        Current_count = 0
        for j in range(len(List_file)):
            if List_file_nonredundant[i] == List_file[j]:
                Current_count += 1
        output_result.write(List_file_nonredundant[i] + "\t" + str(Current_count) + "\n")
