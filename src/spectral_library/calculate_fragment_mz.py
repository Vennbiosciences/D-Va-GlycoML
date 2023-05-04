#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
################################################################################
This script calculates the mass/z for the fragmented ions.
Collection of functions to calculate fragment m/z from the peptide sequence, glycan sequence and ion information.
For peptide backbones, it counts for b and y ions.
For glycan ions, it counts the whole intact peptide plus the mass of glycan.

Created on 18 October 2021.
Modified on 23 February 2022 for five monosaccharides.
Modified on 24 February 2022 for combining amino acid 'I' and 'L' into 'X', to support the fifth monosaccharide.
Modified on 28 February 2022 for handling the fixed modification for amino acid 'C', such as Carbamidomethyl[C].
Modified on 04 June 2022 for calculating precursor ion mass, and annotations.
Modified on 27 August 2022 for handling nine fragmented ions.
################################################################################
"""
__author__ = 'ZLiang'

from pyteomics import mass

# Using special characters to represent five monosaccharides
monosaccharide_component_replacements = {
    "H": "!",
    "N": "@",
    "F": "#",
    "A": "$",
    "G": "%",
}

dumb_reversal = {
    '!': 'B',
    '@': 'O',
    '#': 'U',
    '$': 'I',
    '%': 'L',
}

# Converting special characters back to represent five monosaccharides
character_monosaccharide_replacements = {
    'B': 'H',
    'O': 'N',
    'U': 'F',
    'I': 'A',
    'L': 'G',
}

# Assume the maximum length of peptide is 32, and maximum length of glycan is 18.
MAX_PEPTIDE_LENGTH = 32
MAX_GLYCAN_LENGTH = 18

# Chemical elements for the composition
CHEMICAL_ELEM = 'HCSON'

# Glycan monosaccharides
GLYCAN_MONO = "HNFAG"

# Combine 20 amino acids and 5 monosaccharides together to generate 25 codes, the total number is 26 after plus zero.
# An additional amino acid is 'J', which represents the glycosylated cite for 'N'.
# Ignore 'X' which is used to represent both 'I' and 'L'.
# Should consider the fixed modification for amino acid "C", such as Carbamidomethyl[C].
# The regular expression in mass.py is: "_modX_sequence = re.compile(r'^([^-]+-)?((?:[^A-Z-]*[A-Z])+)(-[^-]+)?$')".
# We Combine 'I' and 'L' into 'X', save the room for the fifth monosaccharide.
# Using 'I' and 'L' to represent monosaccharides, and 'Z' to represent empty one.
# The original one:  amino_acid_codes = "ACDEFGHIJKLMNPQRSTVWY"
amino_acid_codes = "ACDEFGHJKMNPQRSTVWXY"
monosaccharide_codes = "".join(monosaccharide_component_replacements.values())
zero_code = "Z"
amino_acid_monosaccharide_zero_codes = amino_acid_codes + monosaccharide_codes + zero_code
print(amino_acid_monosaccharide_zero_codes)

## Set new dictionary elements for the characters we defined for monosaccharides, 'J', '-'
aa_comp = dict(mass.std_aa_comp)

## go from our new coding back to chemical formula
#   "X": "I" or "L", 113.084064 C(6)H(11)N(1)O(1)S(0)
#   "C": "C" for Carbamidomethyl[C], C(3)H(5)N(1)O(1)S(1) + C(2)H(3)N(1)O(1)
#   "C": "C" for Carbamidomethyl[C], 103.009185 + 57.021464 = 160.030649, C(5)H(8)N(2)O(2)S(1)
#   "H": "!",  Hex H 162.0528234315 20 C(6)O(5)H(10)
#   "N": "@",  HexNAc N 203.079372533 20 C(8)O(5)N(1)H(13)
#   "F": "#",  dHex F 146.0579088094 7 C(6)O(4)H(10)
#   "A": "$",  NeuAc A 291.09541652769997 7 C(11)O(8)N(1)H(17)
#   "G": "%",  NeuGc G 307.09033114979997 7 C(11)O(9)N(1)H(17)
#   "J": "J",  The glycosylated cite for 'N' J 114.042927 C(4)H(6)N(2)O(2)S(0)
#   "Z": "Z",  Zero or empty character 0.0 C(0)

aa_comp["X"] = mass.Composition('C6H11N1O1')  ## chemical formulas for Amino Acid X (I or L)
aa_comp["C"] = mass.Composition('C5H8N2O2S1')  ## chemical formulas for Amino Acid C with Carbamidomethyl[C]
aa_comp["B"] = mass.Composition('C6O5H10')  ## chemical formulas for Hex H
aa_comp["O"] = mass.Composition('C8O5N1H13') ## chemical formulas for HexNAc N
aa_comp["U"] = mass.Composition('C6O4H10')  ## chemical formulas for dHex F
aa_comp["I"] = mass.Composition('C11O8N1H17')  ## chemical formulas for NeuAc A
aa_comp["L"] = mass.Composition('C11O9N1H17')  ## chemical formulas for NeuGc G
aa_comp["J"] = mass.Composition('C4H6N2O2')   ## chemical formulas for J (Glycosylated cite for 'N')
aa_comp["Z"] = mass.Composition('C0')   ## chemical formulas for z (zero or empty)

#  Add mass for the cross ring monosaccharide.
CROSS_RING_MASS = 83.037
# Chemical formulas for cross ring monosaccharide: 'C(4)H(5)N(1)O(1)'
# HCSON = [5, 4, 0, 1, 1]
CROSS_RING_COMP = [5, 4, 0, 1, 1]


class CalculateMZ(object):
    # Calculates the monoisotopic mass of a glycopeptide.
    def calculate_mass(self, glycopeptide_seq, ion_charge):
        """
        :param glycopeptide_seq: string, glycopeptide sequence.
        :param ion_charge: ion charge, integer, such as 1, 2, 3;
        :return: the mass of the precursor ion.
        """
        # Check the input types
        if not isinstance(glycopeptide_seq, (str, )):
            raise ValueError("glycopeptide_seq must be a string")
        if not isinstance(ion_charge, (int,)):
            raise ValueError("ion_charge must be an integer")

        for key, value in dumb_reversal.items():
            glycopeptide_seq = glycopeptide_seq.replace(key, value)

        return mass.calculate_mass(sequence=glycopeptide_seq, charge=ion_charge, aa_comp=aa_comp)


class CalculateFragmentMZ(object):
    # Transfer one hot encodings back into the sequence
    @staticmethod
    def reverse_one_hot_encode(vectors, code):
        """
        read one hot encoding, and convert it into a sequence
        :param vectors: one hot encoding, lists, such as [[0,0,1,...,0],...,[0,1,0,...,0]];
        :param code: rule for the encoding, such as "ACDEFGHJKMNPQRSTVWXY!@#$%Z"
        :return: the corresponding sequence, such as "YKJNSDXSSTRZZZZZZZZZZZZZZZZZZZZZ@@!!!!!!ZZZZZZZZZZ"
        """
        # Check the input types
        if not isinstance(vectors, (list, )):
            raise ValueError("vectors must be a list")
        if not isinstance(code, (str, )):
            raise ValueError("code must be a string")

        letters = []
        for vector in vectors:
            i = vector.index(1)  # get the index of the item which is 1
            letters.append(code[i])
        return "".join(letters)


    # Get fragmented ion mass from one hot encoding
    def get_frag_mz(self, one_hot, ion_position, ion_number, ion_charge):
        """
        :param one_hot: one hot encoding, lists, such as [[0,0,1,...,0],...,[0,1,0,...,0]];
        :param ion_position: integer, ion position, such as 4;
        :param ion_number: ion type of number, b, b$, b-N(1), y, u$, b-N(1), Y0, Y$, Y ions, such as 0, 1, 2, ..., 7, 8;
        :param ion_charge: ion charge, integer, such as 1, 2, 3;
        :return: the mass of the fragment ion, and the total number of different chemical elements.
        """

        # Check the input types
        if not isinstance(one_hot, (list, )):
            raise ValueError("one_hot must be a list")
        if not isinstance(ion_position, (int, )):
            raise ValueError("ion_position must be an integer")
        if not isinstance(ion_number, (int,)):
            raise ValueError("ion_number must be an integer")
        if not isinstance(ion_charge, (int,)):
            raise ValueError("ion_charge must be an integer")

        pep_seq = self.reverse_one_hot_encode(one_hot, amino_acid_monosaccharide_zero_codes)
        #print(pep_seq)
        """
        For b ion, the position is the length of b ion, the length is from left to right;
        For y ion, the position is the "left peptide length" + "length of y ion", the length is from right to left;  
            left_peptide_length = max_peptide_length - len(peptide_name)
            position += left_peptide_length
        For Y ion, the position is the "max_peptide_length" + "length of Y ion", the length is from left to right;
        And ion_type should be 0, 1, 2, ... , 8 which represent b, b$, b-N(1), y, y$, y-N(1), Y, Y0, Y$ ions.
        """
        if ion_number in [0, 1]:
            ion_type = 'b'
            ion_seq = pep_seq[:ion_position]
        elif ion_number == 2:
            # b-N(1), adding "O"
            ion_type = 'b'
            ion_seq = f'{pep_seq[:ion_position]}O'
        elif ion_number in [3, 4]:
            ion_type = 'y'
            # ion_position is the position of y-ion + (MAX_GLYCAN_LENGTH - peptide_length)
            ion_seq = pep_seq[- ion_position - MAX_GLYCAN_LENGTH : MAX_PEPTIDE_LENGTH]
        elif ion_number == 5:
            # y-N(1), adding "O"
            ion_type = 'y'
            ion_seq = f'{pep_seq[- ion_position - MAX_GLYCAN_LENGTH : MAX_PEPTIDE_LENGTH]}O'
        else:
            ion_type = 'y'
            ion_seq = pep_seq[:ion_position]
        #print("ion_position:" + str(ion_position))
        #print("ion_seq:" + ion_seq)

        # Mass of b - ions = Σ(residue masses) + 1(H+)
        # b-ions: m/z = (Σ(residue masses) + 1*z) / z
        # Mass of y - ions = Σ(residue masses) + 19(H2O + H+)
        # y-ions: m / z = (Σ(residue masses) + 1 * z + 18) / z
        for key, value in dumb_reversal.items():
            ion_seq = ion_seq.replace(key, value)
        #print("ion_seq(replace):" + ion_seq)
        mz = mass.calculate_mass(sequence=ion_seq, ion_type=ion_type, charge=int(ion_charge), aa_comp=aa_comp)
        #print(mz)

        chemical_elem_num = [0, 0, 0, 0, 0]
        for aa_mono in ion_seq:
            for key in aa_comp[aa_mono]:
                chemical_elem_num[CHEMICAL_ELEM.find(key)] += aa_comp[aa_mono][key]
        chemical_formulas = "".join(
            f'{CHEMICAL_ELEM[i]}({str(chemical_elem_num[i])})'
            for i in range(len(CHEMICAL_ELEM))
            if chemical_elem_num[i] != 0
        )

        # For the Cases of 1, 4, 7, which contain additional mass annotated as "$"
        # add mass for the cross ring fragment of a HexNAc.
        if ion_number in [1, 4, 7]:
            mz += CROSS_RING_MASS/int(ion_charge)
            # Add composition of HCSON = [5, 4, 0, 1, 1]
            for i in range(len(chemical_elem_num)):
                chemical_elem_num[i] += CROSS_RING_COMP[i]
            chemical_formulas = "".join(
                f'{CHEMICAL_ELEM[i]}({str(chemical_elem_num[i])})'
                for i in range(len(CHEMICAL_ELEM))
                if chemical_elem_num[i] != 0
            )

        return mz, chemical_formulas


    def get_frag_annotation(self, glycopeptide_seq, ion_position, ion_number, ion_charge):
        """
        :param glycopeptide_seq: string, glycopeptide sequence.
        :param ion_position: integer, ion position, such as 4;
        :param ion_number: ion type of number, b, b$, b-N(1), y, u$, b-N(1), Y0, Y$, Y ions, such as 0, 1, 2, ..., 7, 8;
        :param ion_charge: ion charge, integer, such as 1, 2, 3;
        :return: the mass of the fragment ion.
        """

        # Check the input types
        if not isinstance(glycopeptide_seq, (str, )):
            raise ValueError("glycopeptide_seq must be a string")
        if not isinstance(ion_position, (int, )):
            raise ValueError("ion_position must be an integer")
        if not isinstance(ion_number, (int,)):
            raise ValueError("ion_number must be an integer")
        if not isinstance(ion_charge, (int,)):
            raise ValueError("ion_charge must be an integer")

        """
        For b ion, the position is the length of b ion, the length is from left to right;
        For y ion, the position is the "left peptide length" + "length of y ion", the length is from right to left;  
            left_peptide_length = MAX_GLYCAN_LENGTH - len(peptide_name)
            position += left_peptide_length
        For Y ion, the position is the "MAX_GLYCAN_LENGTH" + "length of Y ion", the length is from left to right;
        """
        if ion_number == 0:
            return f'b{str(ion_position)}+{str(ion_charge)}'
        elif ion_number == 1:
            return f'b${str(ion_position)}+{str(ion_charge)}'
        elif ion_number == 2:
            return f'b{str(ion_position)}-N(1)+{str(ion_charge)}'
        elif ion_number in [3, 4, 5]:
            # ion_position is the position of y-ion + (MAX_GLYCAN_LENGTH - peptide_length)
            ion_seq = glycopeptide_seq[- ion_position - MAX_GLYCAN_LENGTH : MAX_PEPTIDE_LENGTH]
            i = 0
            for i in range(len(ion_seq)):
                #print(ion_seq)
                if ion_seq[i] == 'Z':
                    break
            if ion_number == 3:
                return f'y{i}+{str(ion_charge)}'
            elif ion_number == 4:
                return f'y${i}+{str(ion_charge)}'
            else:
                return f'y{i}-N(1)+{str(ion_charge)}'
        elif ion_number == 6:
            return f'Y0+{str(ion_charge)}'
        elif ion_number == 7:
            return f'Y$+{str(ion_charge)}'
        else:
            ion_seq = glycopeptide_seq[MAX_PEPTIDE_LENGTH:ion_position]
            glycan_number = [0, 0, 0, 0, 0]
            for key, value in character_monosaccharide_replacements.items():
                ion_seq = ion_seq.replace(key, value)
            for c in ion_seq:
                glycan_number[GLYCAN_MONO.find(c)] += 1
            glycan_annotation = "".join(
                f'{GLYCAN_MONO[i]}({str(glycan_number[i])})'
                for i in range(len(GLYCAN_MONO))
                if glycan_number[i] != 0
            )
            return f'Y-{glycan_annotation}+{str(ion_charge)}'