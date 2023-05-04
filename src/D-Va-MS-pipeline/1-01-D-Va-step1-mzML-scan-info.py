import math
import os
import sys
from tkinter import filedialog

from pyteomics import mzml
from pyteomics.auxiliary import cvquery


### get file list from a directory
### input (directory)
### output (file list)
def Func_file_list(dir, appendix="ms1"):
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

# input-directory, mzML files
dir = filedialog.askdirectory(initialdir=(os.getcwd()))
os.chdir(dir)
File_list = Func_file_list(dir, "mzML")

# output-info of each scan
with open("mzML-scan-info.txt", "w") as output_result1:
    with open("mzML-scan-info-MS1.txt", "w") as output_result2:
        with open("mzML-scan-info-MS2-HCD45.txt", "w") as output_result3:
            with open("mzML-scan-info-MS2-HCD29.txt", "w") as output_result4:
                with open("mzML-scan-info-MS2-HCD30.txt", "w") as output_result5:
                    output_result1.write("File\tScan\tFile-Scan\tHeader\tRT\n")
                    output_result2.write("File\tScan\tFile-Scan\tHeader\tRT\n")
                    output_result3.write("File\tScan\tFile-Scan\tHeader\tRT\n")
                    output_result4.write("File\tScan\tFile-Scan\tHeader\tRT\n")
                    output_result5.write("File\tScan\tFile-Scan\tHeader\tRT\n")
                    for File_each in File_list:
                        print("processing   ", File_each)
                        for scan in mzml.read(File_each):
                            output_result1.write(File_each.split(".")[0] + "\t")
                            output_result1.write(
                                cvquery(scan, "MS:1000796").split("scan=")[1].replace('"', "")
                                + "\t"
                            )
                            output_result1.write(
                                File_each.split(".")[0]
                                + "-"
                                + cvquery(scan, "MS:1000796").split("scan=")[1].replace('"', "")
                                + "\t"
                            )
                            output_result1.write(cvquery(scan, "MS:1000512") + "\t")
                            output_result1.write(str(cvquery(scan, "MS:1000016")) + "\n")
                            if "p NSI" in cvquery(scan, "MS:1000512"):
                                output_result2.write(File_each.split(".")[0] + "\t")
                                output_result2.write(
                                    cvquery(scan, "MS:1000796").split("scan=")[1].replace('"', "")
                                    + "\t"
                                )
                                output_result2.write(
                                    File_each.split(".")[0]
                                    + "-"
                                    + cvquery(scan, "MS:1000796").split("scan=")[1].replace('"', "")
                                    + "\t"
                                )
                                output_result2.write(cvquery(scan, "MS:1000512") + "\t")
                                output_result2.write(str(cvquery(scan, "MS:1000016")) + "\n")
                            if "hcd45" in cvquery(scan, "MS:1000512"):
                                output_result3.write(File_each.split(".")[0] + "\t")
                                output_result3.write(
                                    cvquery(scan, "MS:1000796").split("scan=")[1].replace('"', "")
                                    + "\t"
                                )
                                output_result3.write(
                                    File_each.split(".")[0]
                                    + "-"
                                    + cvquery(scan, "MS:1000796").split("scan=")[1].replace('"', "")
                                    + "\t"
                                )
                                output_result3.write(cvquery(scan, "MS:1000512") + "\t")
                                output_result3.write(str(cvquery(scan, "MS:1000016")) + "\n")
                            if "hcd29" in cvquery(scan, "MS:1000512"):
                                output_result4.write(File_each.split(".")[0] + "\t")
                                output_result4.write(
                                    cvquery(scan, "MS:1000796").split("scan=")[1].replace('"', "")
                                    + "\t"
                                )
                                output_result4.write(
                                    File_each.split(".")[0]
                                    + "-"
                                    + cvquery(scan, "MS:1000796").split("scan=")[1].replace('"', "")
                                    + "\t"
                                )
                                output_result4.write(cvquery(scan, "MS:1000512") + "\t")
                                output_result4.write(str(cvquery(scan, "MS:1000016")) + "\n")
                            if "hcd30" in cvquery(scan, "MS:1000512"):
                                output_result5.write(File_each.split(".")[0] + "\t")
                                output_result5.write(
                                    cvquery(scan, "MS:1000796").split("scan=")[1].replace('"', "")
                                    + "\t"
                                )
                                output_result5.write(
                                    File_each.split(".")[0]
                                    + "-"
                                    + cvquery(scan, "MS:1000796").split("scan=")[1].replace('"', "")
                                    + "\t"
                                )
                                output_result5.write(cvquery(scan, "MS:1000512") + "\t")
                                output_result5.write(str(cvquery(scan, "MS:1000016")) + "\n")
