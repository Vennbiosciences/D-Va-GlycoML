### Convert MGF format from Rawconvert to pParse

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


## program
## program

# input-directory, mgf files
dir = filedialog.askdirectory(initialdir=(os.getcwd()))
os.chdir(dir)
list_mgf_initial = Func_file_list(dir, "mgf")
list_mgf = []

# ignore formatted MGF files
for i in range(len(list_mgf_initial)):
    if "formatted" not in list_mgf_initial[i]:
        list_mgf.append(list_mgf_initial[i])

# output reformated MGF files from Rawconvert to pParse
for eachfile in list_mgf:
    mgf_filename = eachfile
    mgf_file = open(mgf_filename)
    print("processing now  " + mgf_filename)

    for line in mgf_file:
        with open(mgf_filename.split(".")[0] + "_formatted.mgf", "w") as output:
            list_peak = []
            list_intensity = []
            for line in mgf_file:
                line = line.strip()
                if "BEGIN" in line:
                    list_peak = []
                    list_intensity = []
                elif "TITLE" in line:
                    info_raw = line.split("=")[1].split("\\")[-1].split(".")[0]
                elif "SCANS" in line:
                    info_scan = line.split("=")[1]
                elif "RTINSECONDS" in line:
                    info_RT_sec = "%.4f" % (float(line.split("=")[1]))
                elif "CHARGE" in line:
                    info_Z = line.split("=")[1].split("+")[0]
                elif "PEPMASS" in line:
                    info_MZ = line.split("=")[1]
                elif " " in line and "END" not in line:
                    list_peak.append("%.3f" % (float(line.split(" ")[0])))
                    list_intensity.append("%.1f" % (float(line.split(" ")[1])))
                elif "END" in line:
                    output.write("BEGIN IONS\n")
                    output.write(
                        "TITLE="
                        + info_raw
                        + "."
                        + info_scan
                        + "."
                        + info_scan
                        + "."
                        + info_Z
                        + ".0.dta\n"
                    )
                    output.write("CHARGE=" + info_Z + "+\n")
                    output.write("RTINSECONDS=" + info_RT_sec + "\n")
                    output.write("PEPMASS=" + info_MZ + "\n")
                    for i in range(len(list_peak)):
                        output.write(str(list_peak[i]) + " " + str(list_intensity[i]) + "\n")
                    output.write("2500.000 1.0\n")
                    output.write("END IONS\n")
