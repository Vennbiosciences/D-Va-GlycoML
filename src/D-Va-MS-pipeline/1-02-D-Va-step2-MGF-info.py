### Output MGF infomation, e.g. marker intensity, No of peaks, total intensity

marker_list = [
    138.0550,
    144.0655,
    163.0606,
    168.0661,
    186.0766,
    204.0871,
    366.1400,
    274.0927,
    292.1032,
    243.0174,
    405.0702,
    284.0440,
    446.0968,
]
PPM_fragment = 20

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


### check mass accuracy
### input (mass information)
### output (true or false)
def Func_mass_accuracy(precursor, mass_additional, mass_compare, tolerance, tolerance_unit):
    output = False
    if tolerance_unit == "ppm" and abs(
        float(precursor) - float(mass_additional) - float(mass_compare)
    ) / float(precursor) * 1000000 < float(tolerance):
        output = True
    if tolerance_unit == "Da" and abs(
        float(precursor) - float(mass_additional) - float(mass_compare)
    ) < float(tolerance):
        output = True
    return output


### output base peak
### input (all peaks)
### output (base peak)
def Func_base_peak(list_peak, list_intensity):
    Func_base = 1
    Func_intensity_top = 1
    for i in range(len(list_peak)):
        if list_intensity[i] > Func_intensity_top:
            Func_intensity_top = list_intensity[i]
            Func_base = "%.4f" % (list_peak[i])
    return Func_base, Func_intensity_top


### output targeted marker intensity
### input (all peaks, marker info)
### output (marker intensity within PPM)
def Func_marker_intensity(list_peak, list_intensity, marker, PPM):
    Func_intensity = 0
    for i in range(len(list_peak)):
        if (
            Func_mass_accuracy(list_peak[i], 0, marker, PPM, "ppm")
            and list_intensity[i] > Func_intensity
        ):
            Func_intensity = list_intensity[i]
    return Func_intensity


## program
## program

# input-directory, mgf files
dir = filedialog.askdirectory(initialdir=(os.getcwd()))
os.chdir(dir)
list_mgf = Func_file_list(dir, "mgf")

marker_base_range = 10

# output 1)info of each spectrum; 2)summary of each MGF
with open("MGF-info.txt", "w") as output:
    with open("MS2-TIC-204.txt", "w") as output2:
        output.write(
            "RAW\tSpectra\tScan\tRT_min\tRT_sec\tMZ\tM\tZ\tNo_Peak\tIntensity_total\tBase_Peak\tBase_Peak_Int\t138area_Base_Int\t138area_Base%\t204area_Base_Int\t204area_Base%"
        )
        output2.write(
            "RAW\tMS2_TIC_total\tMS2_TIC_138\tMS2_TIC_138&204\tMS2_TIC_163\tMS2_TIC_274\t138%\t138&204%\t163%\t274%\n"
        )
        for i in range(len(marker_list)):
            output.write("\tInt_" + str(marker_list[i]))
        output.write("\n")
        for eachfile in list_mgf:
            mgf_filename = eachfile
            mgf_file = open(mgf_filename)
            print("processing now  " + mgf_filename)
            MS2_TIC_total = 0
            MS2_TIC_138 = 0
            MS2_TIC_138and204 = 0
            MS2_TIC_163 = 0
            MS2_TIC_274 = 0
            for line in mgf_file:
                line = line.strip()
                if "BEGIN" in line:
                    list_peak = []
                    list_intensity = []
                    list_peak_138area = []
                    list_intensity_138area = []
                    list_peak_204area = []
                    list_intensity_204area = []
                    list_peak_138area.append(0)
                    list_intensity_138area.append(0)
                    list_peak_204area.append(0)
                    list_intensity_204area.append(0)
                elif "TITLE" in line:
                    info_spectra = line.split(" ")[0].split("=")[1]
                    info_raw = line.split(" ")[0].split("=")[1].split(".")[0]
                    info_scan = line.split(" ")[0].split("=")[1].split(".")[1]
                elif "RTINSECONDS" in line:
                    info_RT_sec = "%.4f" % (float(line.split("=")[1]))
                    info_RT_min = "%.4f" % (float(info_RT_sec) / 60)
                elif "PEPMASS" in line:
                    info_MZ = "%.4f" % (float(line.split(" ")[0].split("=")[1]))
                elif "CHARGE" in line:
                    info_Z = "%.4f" % (float(line.split("=")[1].split("+")[0]))
                elif "END" not in line:
                    info_M = "%.4f" % (float(info_MZ) * float(info_Z) - float(info_Z) * 1.0073)
                    list_peak.append(float(line.split(" ")[0]))
                    list_intensity.append(float(line.split(" ")[1]))
                    if (
                        marker_list[0] - marker_base_range
                        <= (float(line.split(" ")[0]))
                        <= marker_list[0] + marker_base_range
                    ):
                        list_peak_138area.append(float(line.split(" ")[0]))
                        list_intensity_138area.append(float(line.split(" ")[1]))
                    if (
                        marker_list[5] - marker_base_range
                        <= (float(line.split(" ")[0]))
                        <= marker_list[5] + marker_base_range
                    ):
                        list_peak_204area.append(float(line.split(" ")[0]))
                        list_intensity_204area.append(float(line.split(" ")[1]))
                elif "END" in line:
                    output.write(
                        info_raw
                        + "\t"
                        + info_spectra
                        + "\t"
                        + info_scan
                        + "\t"
                        + info_RT_min
                        + "\t"
                        + info_RT_sec
                        + "\t"
                        + info_MZ
                        + "\t"
                        + info_M
                        + "\t"
                        + info_Z
                        + "\t"
                        + str(len(list_peak))
                        + "\t"
                        + str(sum(list_intensity))
                        + "\t"
                        + str(Func_base_peak(list_peak, list_intensity)[0])
                        + "\t"
                        + str(Func_base_peak(list_peak, list_intensity)[1])
                        + "\t"
                        + str(Func_base_peak(list_peak_138area, list_intensity_138area)[1])
                        + "\t"
                        + str(
                            Func_marker_intensity(list_peak, list_intensity, 138.0550, PPM_fragment)
                            / Func_base_peak(list_peak_138area, list_intensity_138area)[1]
                        )
                        + "\t"
                        + str(Func_base_peak(list_peak_204area, list_intensity_204area)[1])
                        + "\t"
                        + str(
                            Func_marker_intensity(list_peak, list_intensity, 204.0871, PPM_fragment)
                            / Func_base_peak(list_peak_204area, list_intensity_204area)[1]
                        )
                    )
                    for i in range(len(marker_list)):
                        output.write(
                            "\t"
                            + str(
                                Func_marker_intensity(
                                    list_peak, list_intensity, marker_list[i], PPM_fragment
                                )
                            )
                        )
                    output.write("\n")
                    MS2_TIC_total += sum(list_intensity)
                    if (
                        Func_marker_intensity(
                            list_peak, list_intensity, marker_list[0], PPM_fragment
                        )
                        / Func_base_peak(list_peak_138area, list_intensity_138area)[1]
                        + Func_marker_intensity(
                            list_peak, list_intensity, marker_list[5], PPM_fragment
                        )
                        / Func_base_peak(list_peak_204area, list_intensity_204area)[1]
                    ) > 1.5:
                        MS2_TIC_138and204 += sum(list_intensity)
                    if (
                        Func_marker_intensity(
                            list_peak, list_intensity, marker_list[0], PPM_fragment
                        )
                        / Func_base_peak(list_peak_138area, list_intensity_138area)[1]
                        > 0.9
                    ):
                        MS2_TIC_138 += sum(list_intensity)
                        if (
                            Func_marker_intensity(
                                list_peak, list_intensity, marker_list[2], PPM_fragment
                            )
                            / (
                                Func_marker_intensity(
                                    list_peak, list_intensity, marker_list[5], PPM_fragment
                                )
                                + 1
                            )
                            > 0.1
                        ):
                            MS2_TIC_163 += sum(list_intensity)
                        if (
                            Func_marker_intensity(
                                list_peak, list_intensity, marker_list[7], PPM_fragment
                            )
                            / (
                                Func_marker_intensity(
                                    list_peak, list_intensity, marker_list[5], PPM_fragment
                                )
                                + 1
                            )
                            > 0.1
                        ):
                            MS2_TIC_274 += sum(list_intensity)
            output2.write(
                mgf_filename
                + "\t"
                + str(MS2_TIC_total)
                + "\t"
                + str(MS2_TIC_138)
                + "\t"
                + str(MS2_TIC_138and204)
                + "\t"
                + str(MS2_TIC_163)
                + "\t"
                + str(MS2_TIC_274)
                + "\t"
                + str("%.3f" % (float(MS2_TIC_138 / MS2_TIC_total)))
                + "\t"
                + str("%.3f" % (float(MS2_TIC_138and204 / MS2_TIC_total)))
                + "\t"
                + str("%.3f" % (float(MS2_TIC_163 / MS2_TIC_total)))
                + "\t"
                + str("%.3f" % (float(MS2_TIC_274 / MS2_TIC_total)))
                + "\n"
            )
