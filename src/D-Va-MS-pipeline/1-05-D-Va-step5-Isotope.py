import os
from tkinter import filedialog

from brainpy import isotopic_variants


### peptide composition
### input (peptide sequence)
### output (peptide composition)
def Func_Seq_Pep(Peptide):
    Output_Result = [0, 0, 0, 0, 0]
    for i in range(len(Composition_AA_List)):
        Count_AA = Peptide.count(Composition_AA_List[i])
        for j in range(len(Output_Result)):
            Output_Result[j] += Composition_AA[i][j] * Count_AA
    Output_Result[1] += 2
    Output_Result[2] += 1
    return Output_Result


### glycan composition
### input (glycan sequence)
### output (glycan composition)
def Func_Seq_Gly(Glycan):
    Output_Result = [0, 0, 0, 0, 0]
    for i in range(len(Composition_Glycan_List)):
        if f"{Composition_Glycan_List[i]}(" in Glycan:
            Count_Glycan = Glycan.split(f"{Composition_Glycan_List[i]}(")[1].split(")")[0]
            for j in range(len(Output_Result)):
                Output_Result[j] += Composition_Glycan[i][j] * int(Count_Glycan)
    return Output_Result


### composition for Skyline
### input (composition in number)
### output (composition for Skyline)
def Func_format_composition(Composition):
    List_element = ["C", "H", "O", "N", "S"]
    Output_Result = ""
    for i in range(len(Composition)):
        Output_Result += List_element[i]
        Output_Result += str(Composition[i])
    return Output_Result


### calculate MW
### input (composition)
### output (NW)
def Func_calculate_MW(Composition):
    MW_element = [12.000000, 1.007825, 15.994914, 14.003073, 31.972070]
    return sum(Composition[i] * MW_element[i] for i in range(len(Composition)))


### export list of column in txt file
### input (txt file, column header)
### output (column data)
def Func_read_txt_column(Txt_filename, Column_name):
    Txt_file = [line.strip() for line in open(Txt_filename).read().splitlines()]
    Txt_data = [line.split("\t") for line in Txt_file]
    for i in range(len(Txt_data[0])):
        if Txt_data[0][i] == Column_name:
            Column_target = i
    return [Txt_data[i][Column_target] for i in range(1, len(Txt_data))]


### calculate isotope M0
### input (composition)
### output (M0%)
def Func_isotope_M0(Element_list):
    glycopeptide = {
        "C": Element_list[0],
        "H": Element_list[1],
        "O": Element_list[2],
        "N": Element_list[3],
        "S": Element_list[4],
    }
    theoretical_isotopic_cluster = isotopic_variants(glycopeptide, npeaks=20, charge=1)
    return float(theoretical_isotopic_cluster[0].intensity)


### calculate isotope M1
### input (composition)
### output (M1%)
def Func_isotope_M1(Element_list):
    glycopeptide = {
        "C": Element_list[0],
        "H": Element_list[1],
        "O": Element_list[2],
        "N": Element_list[3],
        "S": Element_list[4],
    }
    theoretical_isotopic_cluster = isotopic_variants(glycopeptide, npeaks=20, charge=1)
    return float(theoretical_isotopic_cluster[1].intensity)


## program
## program

# Composition of each AA (A to Z, 26 total) (Element number of C, H, O, N, S)
Composition_AA = []
for i in range(26):
    Composition_AA.append([])
Composition_AA_List = [
    "A",
    "B",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "J",
    "K",
    "L",
    "M",
    "N",
    "O",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "U",
    "V",
    "W",
    "X",
    "Y",
    "Z",
]
Composition_AA[0] = [3, 5, 1, 1, 0]  # A
Composition_AA[1] = [0, 0, 0, 0, 0]  # B
Composition_AA[2] = [5, 8, 2, 2, 1]  # C with Carbamidomethyl
Composition_AA[3] = [4, 5, 3, 1, 0]  # D
Composition_AA[4] = [5, 7, 3, 1, 0]  # E
Composition_AA[5] = [9, 9, 1, 1, 0]  # F
Composition_AA[6] = [2, 3, 1, 1, 0]  # G
Composition_AA[7] = [6, 7, 1, 3, 0]  # H
Composition_AA[8] = [6, 11, 1, 1, 0]  # I
Composition_AA[9] = [4, 6, 2, 2, 0]  # J
Composition_AA[10] = [6, 12, 1, 2, 0]  # K
Composition_AA[11] = [6, 11, 1, 1, 0]  # L
Composition_AA[12] = [5, 9, 1, 1, 1]  # M
Composition_AA[13] = [4, 6, 2, 2, 0]  # N
Composition_AA[14] = [0, 0, 0, 0, 0]  # O
Composition_AA[15] = [5, 7, 1, 1, 0]  # P
Composition_AA[16] = [5, 8, 2, 2, 0]  # Q
Composition_AA[17] = [6, 12, 1, 4, 0]  # R
Composition_AA[18] = [3, 5, 2, 1, 0]  # S
Composition_AA[19] = [4, 7, 2, 1, 0]  # T
Composition_AA[20] = [0, 0, 0, 0, 0]  # U
Composition_AA[21] = [5, 9, 1, 1, 0]  # V
Composition_AA[22] = [11, 10, 1, 2, 0]  # W
Composition_AA[23] = [6, 11, 1, 1, 0]  # X
Composition_AA[24] = [9, 9, 2, 1, 0]  # Y
Composition_AA[25] = [0, 0, 0, 0, 0]  # Z

# Composition of each glycan (Element number of C, H, O, N, S)
Composition_Glycan = []
for i in range(9):
    Composition_Glycan.append([])
Composition_Glycan_List = ["H", "N", "F", "A", "G", "KDN", "X", "HA", "aH"]
Composition_Glycan[0] = [6, 10, 5, 0, 0]  # Hex
Composition_Glycan[1] = [8, 13, 5, 1, 0]  # HexNAc
Composition_Glycan[2] = [6, 10, 4, 0, 0]  # Fuc
Composition_Glycan[3] = [11, 17, 8, 1, 0]  # NeuAc
Composition_Glycan[4] = [11, 17, 9, 1, 0]  # NeuGc
Composition_Glycan[5] = [9, 14, 8, 0, 0]  # KDN
Composition_Glycan[6] = [5, 8, 4, 0, 0]  # Xyl
Composition_Glycan[7] = [6, 8, 6, 0, 0]  # HexA
Composition_Glycan[8] = [6, 13, 5, 1, 0]  # aH

# input-txt files
Input_filename = filedialog.askopenfilename(
    initialdir=(os.getcwd()), filetypes=[("input data", ".txt")], title=("input data")
)
print(Input_filename)

Spectra_name = Func_read_txt_column(Input_filename, "GlySpec")
Sequence_peptide = Func_read_txt_column(Input_filename, "Peptide")
Sequence_glycan = Func_read_txt_column(Input_filename, "GlycanComposition")
TIC_mono = Func_read_txt_column(Input_filename, "MonoArea")
TIC_all = Func_read_txt_column(Input_filename, "IsotopeArea")

Element_peptide = []
Element_glycan = []
Element_glypep = []

# output composition of each glycopeptide
for i in range(len(Sequence_peptide)):
    Element_peptide_current = Func_Seq_Pep(Sequence_peptide[i])
    Element_glycan_current = Func_Seq_Gly(Sequence_glycan[i])
    Element_glypep_current = [0, 0, 0, 0, 0]
    for j in range(len(Element_peptide_current)):
        Element_glypep_current[j] = Element_peptide_current[j] + Element_glycan_current[j]
    Element_peptide.append(Element_peptide_current)
    Element_glycan.append(Element_glycan_current)
    Element_glypep.append(Element_glypep_current)

# check the isotope pass threshold for each glycopeptide
with open(Input_filename.split(".")[0] + "-step5-output.txt", "w") as Output:
    Output.write(
        "GlySpec\tPeptide\tGlycan\tGlyPep_Composition\tGlyPep_MW\tPep_Composition\tPep_MW\tGly_Composition\tGly_MW\tMonoArea\tIsotopeArea\tM0%\tM1%\tD-Va-Step5-check\n"
    )
    for i in range(len(Sequence_peptide)):
        current_mono_check = 1
        if (Func_isotope_M1(Element_glypep[i]) - Func_isotope_M0(Element_glypep[i])) > 0.05:
            if (float(TIC_mono[i]) / (float(TIC_all[i]) + 1)) > (
                (Func_isotope_M1(Element_glypep[i]) + (Func_isotope_M0(Element_glypep[i]))) / 2
            ):
                current_mono_check = 0
        Output.write(
            Spectra_name[i]
            + "\t"
            + Sequence_peptide[i]
            + "\t"
            + Sequence_glycan[i]
            + "\t"
            + str(Func_format_composition(Element_glypep[i]))
            + "\t"
            + str(Func_calculate_MW(Element_glypep[i]))
            + "\t"
            + str(Func_format_composition(Element_peptide[i]))
            + "\t"
            + str(Func_calculate_MW(Element_peptide[i]))
            + "\t"
            + str(Func_format_composition(Element_glycan[i]))
            + "\t"
            + str(Func_calculate_MW(Element_glycan[i]))
            + "\t"
            + str(TIC_mono[i])
            + "\t"
            + str(TIC_all[i])
            + "\t"
            + str(Func_isotope_M0(Element_glypep[i]))
            + "\t"
            + str(Func_isotope_M1(Element_glypep[i]))
            + "\t"
            + str(current_mono_check)
            + "\n"
        )

Result_check = Func_read_txt_column(
    Input_filename.split(".")[0] + "-step5-output.txt", "D-Va-Step5-check"
)
Result_check.insert(0, "1")
# output psm pass threshold for each glycopeptide
with open(Input_filename.split(".")[0] + "-step5-psm.txt", "w") as Output:
    Extract_file = [line.strip() for line in open(Input_filename).read().splitlines()]
    for i in range(len(Extract_file)):
        if Result_check[i] == "1":
            Output.write(Extract_file[i])
            Output.write("\n")
