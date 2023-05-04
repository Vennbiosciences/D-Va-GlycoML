import math
import os
import sys
from tkinter import filedialog


### get file list from a directory
### input (directory)
### ouput (file list)
def Func_file_list(dir, appendix="mgf"):
    import os

    file_num = 0
    output_list = []
    file_all = os.listdir(dir)
    for line in file_all:
        filepath = os.path.join(dir, line)
        if filepath.split(".")[-1] == appendix:
            output_list.append(filepath.split("\\")[-1])
            print(filepath.split("\\")[-1])
    return output_list


### check current scan is HCD29
### input (all scan & current scan)
### output (true or false)
def Func_check_HCD29(info_raw, info_scan, list_raw, list_scan):
    return any(info_raw == list_raw[i] and info_scan == list_scan[i] for i in range(len(list_raw)))


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

# input-directory, mgf files
dir = filedialog.askdirectory(initialdir=(os.getcwd()))
os.chdir(dir)
list_mgf = Func_file_list(dir, "mgf")

mzML_scan = filedialog.askopenfilename(
    initialdir=(os.getcwd()), filetypes=[("input data", ".txt")], title=("input data")
)
print(mzML_scan)
list_raw = Func_read_txt_column(mzML_scan, "File")
list_scan = Func_read_txt_column(mzML_scan, "Scan")

# output new mgf file with targetd scan
for eachfile in list_mgf:
    mgf_filename = eachfile
    mgf_file = open(mgf_filename)

    with open(mgf_filename.split(".")[0] + "_selected.mgf", "w") as output:
        list_peak = []
        list_intensity = []
        list_raw_current = []
        list_scan_current = []
        for i in range(len(list_raw)):
            if (
                mgf_filename.split(".")[0].replace("_HCDFT", "") in list_raw[i]
                or mgf_filename.split(".")[0].replace("_formatted", "") in list_raw[i]
            ):
                list_raw_current.append(list_raw[i])
                list_scan_current.append(list_scan[i])
        print("processing   ", eachfile)
        for line in mgf_file:
            line = line.strip()
            if "BEGIN" in line:
                list_peak = []
                list_intensity = []
            elif "TITLE" in line:
                info_title = line.split(" ")[0].split("=")[1]
                info_raw = line.split(" ")[0].split("=")[1].split(".")[0]
                info_scan = line.split(" ")[0].split("=")[1].split(".")[1]
            elif "RTINSECONDS" in line:
                info_RT_sec = "%.4f" % (float(line.split("=")[1]))
            elif "CHARGE" in line:
                info_Z = line.split("=")[1].split("+")[0]
            elif "PEPMASS" in line:
                info_MZ = line.split(" ")[0].split("=")[1]
            elif " " in line and "END" not in line:
                list_peak.append("%.3f" % (float(line.split(" ")[0])))
                list_intensity.append("%.1f" % (float(line.split(" ")[1])))
            elif "END" in line:
                if Func_check_HCD29(info_raw, info_scan, list_raw_current, list_scan_current):
                    output.write("BEGIN IONS\n")
                    output.write("TITLE=" + info_title + "\n")
                    output.write("CHARGE=" + info_Z + "+\n")
                    output.write("RTINSECONDS=" + info_RT_sec + "\n")
                    output.write("PEPMASS=" + info_MZ + "\n")
                    for i in range(len(list_peak)):
                        output.write(str(list_peak[i]) + " " + str(list_intensity[i]) + "\n")
                    output.write("2500.000 1.0\n")
                    output.write("END IONS\n")
