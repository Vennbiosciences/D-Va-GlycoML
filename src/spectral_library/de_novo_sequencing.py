#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
#####################################################################################################
This script generates a linear array for the glycan from an input spectrum.
By using glycan de novo sequencing, a set of linearized glycan will output as candidates.

Created on 1 October 2021.
Modified on 18 February 2022 for five monosaccharides.
Modified on 21 February 2022 for double-digit monosaccharides.
#####################################################################################################
"""
__author__ = 'ZLiang'


monosaccharide_number_replacement = {
    "0100000000": "H",
    "0001000000": "N",
    "0000010000": "F",
    "0000000100": "A",
    "0000000001": "G"
}


class DeNovoSequencing(object):

    @staticmethod
    def count_total(monosaccharide_numbers):
        """
        Summarize the numbers for the monosaccharides.
        Should consider double-digit monosaccharides, here we count every two digits as a number.
        :param monosaccharide_numbers: the sequence of numbers for each monosaccharide, such as '081204010000'
        :return: total number of the monosaccharides. Such as 8+12+4+1=25.
        Note: Should not be 8+1+2+4+1 = 16.
        """

        # Check the input types
        if not isinstance(monosaccharide_numbers, (str, )):
            raise ValueError("monosaccharide_numbers must be a string")

        mono_length = len(monosaccharide_numbers)
        assert (mono_length % 2 == 0), "Wrong Input Monosaccharides for count_total!"
        if monosaccharide_numbers.isdigit():
            # transfer an integer into a list
            numbers = list(monosaccharide_numbers)
        else:
            return "The sequence is not a number!"

        return sum(int(numbers[i] + numbers[i+1]) for i in range(0, mono_length, 2))


    @staticmethod
    def glycan_list_to_str(glycan_list):
        """
        Parse a glycan list such as "HHNFN" into a string "22100", which represents the number for "HNFA".
        Also, should consider double-digit monosaccharides.
        :param glycan_list: A string of list with [character]+, such as "HHNFN"
        :return: A string like "0202010000"
        """

        # Check the input types
        if not isinstance(glycan_list, (str, )):
            raise ValueError("glycan_list must be a string")

        glycan_mono = "HNFAG"
        glycan_number = ['0', '0', '0', '0', '0']

        for i in range(len(glycan_list)):
            glycan_number[glycan_mono.find(glycan_list[i])] = str(int(glycan_number[glycan_mono.find(glycan_list[i])]) + 1)
        # should add 0 to the front of glycan number if it is only one digit
        for i in range(len(glycan_number)):
            if len(glycan_number[i]) == 1:
                glycan_number[i] = '0' + glycan_number[i]

        return "".join(glycan_number)


    @staticmethod
    def glycan_contain(glycan_former, glycan_current):
        """
        Should consider double-digit monosaccharides.
        Judge whether a current glycan list such as "1002010000" contains a former glycan list such as "0902010000"
        :param glycan_former: A former glycan list, such as "0902010000"
        :param glycan_former: A current glycan list, such as "1002010000"
        :return: True or False
        """

        # Check the input types
        if not isinstance(glycan_former, (str, )):
            raise ValueError("glycan_former must be a string")
        if not isinstance(glycan_current, (str, )):
            raise ValueError("glycan_current must be a string")

        assert len(glycan_former) == len(glycan_current)
        glycan_length = len(glycan_former)
        assert (glycan_length % 2 == 0), "Wrong Input Monosaccharides for glycan_contain!"
        return all(
            int (glycan_former[i] + glycan_former[i+1]) <= int(glycan_current[i] + glycan_current[i+1])
            for i in range(0, glycan_length, 2)
        )


    @staticmethod
    def glycan_minus(glycan_former, glycan_current):
        """
        Should call glycan_contain() before this function, to make sure glycan_current contains glycan_former.
        Should consider double-digit monosaccharides, and glycan current should contains glycan former.
        Use a current glycan list such as "1002010000" minuses a former glycan list such as "0902010000"
        :param glycan_former: A former glycan list, such as "0902010000"
        :param glycan_former: A current glycan list, such as "1002010000"
        :return: A string like "0100000000"
        """

        # Check the input types
        if not isinstance(glycan_former, (str, )):
            raise ValueError("glycan_former must be a string")
        if not isinstance(glycan_current, (str, )):
            raise ValueError("glycan_current must be a string")

        assert len(glycan_former) == len(glycan_current)
        glycan_length = len(glycan_former)
        assert (glycan_length % 2 == 0), "Wrong Input Monosaccharides for glycan_minus!"
        # Minus all the letters
        glycan_after_minus = list(range(glycan_length // 2))
        # i starts from 0 to glycan_length/2, j starts from 0 to glycan_length
        for i, j in enumerate(range(0, glycan_length, 2)):
            glycan_after_minus[i] = str(int(glycan_current[j] + glycan_current[j+1]) - int(glycan_former[j] + glycan_former[j+1]))
        # should add 0 to the front of glycan number if it is only one digit
        for i in range(len(glycan_after_minus)):
            if len(glycan_after_minus[i]) == 1:
                glycan_after_minus[i] = '0' + glycan_after_minus[i]

        return "".join(glycan_after_minus)


    @staticmethod
    def glycan_str_to_list(glycan_str):
        """
        Should consider double-digit monosaccharides, and only return a list from the a set of combinatorial lists.
        Parse a glycan string such as "0002010000" into a set of combinatorial lists, such as {"NNA", "NAN", "ANN"}, and
        only return the list "NNA"
        :param glycan_str: A string of glycan numbers [digit]+, such as "0002010000",
        :return: A list like "NNA".
        """

        # Check the input types
        if not isinstance(glycan_str, (str, )):
            raise ValueError("glycan_former must be a string")

        glycan_list = []
        glycan_mono = "HNFAG"
        glycan_length = len(glycan_str)
        assert (glycan_length % 2 == 0), "Wrong Input Monosaccharides for glycan_str_to_list!"
        # Should consider double-digit monosaccharides, which should visit glycan_mono[i//2]
        for i in range(0, glycan_length, 2):
            mono_number = int(glycan_str[i] + glycan_str[i+1]) - 0
            for _ in range(mono_number):
                glycan_list.append(glycan_mono[i//2])

        # Should we consider the combinations of the glycan?
        #glycan_set_list = permutations(glycan_list)
        #glycan_no_repeat_set = list(set(glycan_set_list))
        #return glycan_no_repeat_set
        return "".join(glycan_list)


    def de_novo(self, glycan_precursor, glycan_intensity_dic):
        """
        Generate a set of linearized glycan sequences, from glycan_intensity_dictionary.
        Here we use 2D lists to store the results of de novo sequencing.
        Time complexity should be O(n^3).
        For example:
            PEPTIDE=LSSJSTKK
            Glycan=HHHHHNNNNGG
            Glycan String=0504000002
        :param glycan_precursor: the number representation of glycan precursor, such as "0504000002"
        :param glycan_intensity_dic: a dictionary contains composition and intensity
            Such as {'0001000000': 178721.8, '0002000000': 24745.4, '0102000000': 1714.1, '0404000000': 23748.3,
                '0202000000': 46100.9, '0504000000': 117434.6, '0404000001': 18476.5, '0302000000': 40545.4,
                '0203000000': 3935.4, '0504000001': 25856.8, '0303000000': 76999.1, '0403000000': 233672.0,
                '0303000001': 74035.4, '0403000001': 312205.2}
        :return: glycan de novo string tuple lists, such as [('NNHHNHHANHNHAA', 954586.3), ('NNHHHNHANHNHAA', 954473.8)]
        """

        # Check the input types
        if not isinstance(glycan_precursor, (str, )):
            raise ValueError("glycan_precursor must be a string")
        if not isinstance(glycan_intensity_dic, (dict, )):
            raise ValueError("glycan_precursor must be a dictionary")

        # First step, count total number of glycan, then generate a list with the length of it.
        # Here we use additional one to store the whole sequence
        glycan_precursor_length = self.count_total(glycan_precursor) + 1
        # Total monosaccharides, get the maximum ions sequence number, index [0] will be an empty list here.
        glycan_lists = [[] for _ in range(glycan_precursor_length)]

        # Second step, count the total number of monosaccharides for the ions
        for key, value in glycan_intensity_dic.items():
            ion_number = self.count_total(key)
            ion_dic = {key: value}
            glycan_lists[ion_number].append(ion_dic)

        for i in range(len(glycan_lists)):
            print(glycan_lists[i])

        # Third step, generate glycan linear ladders based on the former ones
        # Should add a start ladder, which might not be the glycan_lists[1], in this case, should generate the sequence.
        glycan_ladders = [[] for _ in range(glycan_precursor_length)]
        first_ladder = True
        # Base case for the first ladder 0, add {"0000000000": 0.0} in it, as an initialization.
        glycan_ladders[0].append({"": 0.0})
        # Base case for the last glycan list:
        # If there is no glycan list in the last one, add {"glycan_precursor": 0.0} in it, as an end;
        # Otherwise, do nothing.
        if len(glycan_lists[glycan_precursor_length - 1]) == 0:
            glycan_lists[glycan_precursor_length - 1].append({glycan_precursor: 0.0})
        #print(glycan_lists)
        for i in range(1, glycan_precursor_length):
            # Base case for the first ladder, might contain several monosaccharides.
            if first_ladder:
                for j in range(len(glycan_lists[i])):
                # Base case for the first monosaccharide
                    for key, value in glycan_lists[i][j].items():
                        glycan_ladder_list = self.glycan_str_to_list(key)
                        new_key = "".join(glycan_ladder_list)
                        glycan_ladders[i].append({new_key: value})
                # Should pad the former ladders, but we ignore the process here, since the final step will contain it.
                first_ladder = False
            # Dynamic programming for generating ladders by top-down methods, from Y-ions to B-ions
            else:
                for j in range(len(glycan_lists[i])):
                    # Back tracing the former ladder one by one until find a successful prefix
                    # Check each element in the former ladder
                    for key, value in glycan_lists[i][j].items():
                        former_ladder = i - 1
                        find = False
                        while not find and former_ladder >= 0:
                            for k in range(len(glycan_ladders[former_ladder])):
                                for former_key, former_value in glycan_ladders[former_ladder][k].items():
                                    former_glycan = self.glycan_list_to_str(former_key)
                                    if self.glycan_contain(former_glycan, key):
                                        find = True
                                        glycan_remains = self.glycan_minus(former_glycan, key)
                                        new_key = former_key + self.glycan_str_to_list(glycan_remains)
                                        new_value = round(former_value + value, 1)
                                        glycan_ladders[i].append({new_key: new_value})
                            former_ladder -= 1

        for i in range(len(glycan_ladders)):
            print(glycan_ladders[i])

        # Use a set to store the last glycan_ladders, then sorted by intensity
        glycan_de_novo = {}
        for i in range(len(glycan_ladders[glycan_precursor_length - 1])):
            for key, value in glycan_ladders[glycan_precursor_length - 1][i].items():
                glycan_de_novo[key] = value
        #for g in range(glycan_precursor_length):
        #    print(glycan_ladders[g])
        return sorted(
            glycan_de_novo.items(), key=lambda kv: (kv[1], kv[0]), reverse=True
        )


    def de_novo_gap(self, glycan, glycan_intensity_dic):
        """
        Generate a set of linearized glycan sequences, from glycan_intensity_dictionary.
        We use 2D lists to store the results of de novo sequencing.
        Time complexity should be O(n^3).
        Here we consider the gap between the ladders:
            if the gap is 0, we do not continue to search the former ladder,
            otherwise, we need to search the former ladders plus the gap,
                for example, current_ladder = 4, gap = 1, we need to search ladder = 3 and ladder = 4
        :param glycan: the number representation of glycans, such as "5410"
        :param glycan_intensity_dic: a dictionary contains composition and intensity
            Such as {'2410': 400.3, '1300': 589.3, '2310': 491.9, '4300': 413.3, '0100': 3591.5,
            '4310': 465.1, '1200': 648.3, '2200': 1793.2, '3200': 1706.6}
        :return: original amino acid name string
        """

        # Check the input types
        if not isinstance(glycan, (str, )):
            raise ValueError("glycan must be a string")
        if not isinstance(glycan_intensity_dic, (dict, )):
            raise ValueError("glycan_precursor must be a dictionary")

        # First step, count total number of glycan, and generate a list with the length of it.
        glycan_length = self.count_total(glycan)
        # total monosaccharides, get the maximum ions sequence number, index [0] will be an empty list here.
        glycan_lists = [[] for _ in range(glycan_length)]

        # Second step, count the total number of monosaccharides for the ions
        for key, value in glycan_intensity_dic.items():
            ion_number = self.count_total(key)
            ion_dic = {key: value}
            glycan_lists[ion_number].append(ion_dic)

        # Third step, generate glycan linear ladders based on the former ones
        glycan_ladders = [[] for _ in range(glycan_length)]
        glycan_gaps = [[] for _ in range(glycan_length)]
        for i in range(1, glycan_length):
            for j in range(len(glycan_lists[i])):
                # Base case for only one monosaccharide
                if i == 1:
                    for key, value in glycan_lists[i][j].items():
                        monosaccharide = monosaccharide_number_replacement.get(key)
                        glycan_ladders[i].append({monosaccharide: value})
                        glycan_gaps[i].append(0)
                # Dynamic programming for generating ladders by top-down methods, from Y-ions to B-ions
                else:
                    # Back tracing the former ladder one by one until find a successful prefix
                    # Check each elements in the former ladder
                    for key, value in glycan_lists[i][j].items():
                        former_ladder = i - 1
                        gap = 0
                        find = False
                        contain = 0
                        while (not find or contain > 0) and former_ladder > 0:
                            # Assume there is no gap, set contain to 0
                            contain = 0
                            for k in range(len(glycan_ladders[former_ladder])):
                                for former_key, former_value in glycan_ladders[former_ladder][k].items():
                                    former_glycan = self.glycan_list_to_str(former_key)
                                    if self.glycan_contain(former_glycan, key):
                                        find = True
                                        if glycan_gaps[former_ladder][k] > contain:
                                            contain = glycan_gaps[former_ladder][k]
                                        glycan_remains = self.glycan_minus(former_glycan, key)
                                        glycan_ladder_list = self.glycan_str_to_list(glycan_remains)
                                        new_key = former_key + "".join(glycan_ladder_list)
                                        new_value = round(former_value + value, 1)
                                        glycan_ladders[i].append({new_key: new_value})
                                        glycan_gaps[i].append(gap)
                                        '''
                                        for glycan_ladder_list in glycan_ladder_lists:
                                            print(former_key)
                                            new_key = former_key + "".join(glycan_ladder_list)
                                            print(new_key)
                                            print(value)
                                            new_value = round(former_value + value, 1)
                                            glycan_ ladders[i].append({new_key: new_value})
                                            glycan_gaps[i].append(gap)
                                            print(glycan_ladders[i])
                                        '''
                            former_ladder -= 1
                            gap += 1

        for g in range(glycan_length):
            print(glycan_ladders[g])
        for p in range(glycan_length):
            print(glycan_gaps[p])
        return glycan_lists


    def de_novo_dic(self, glycan, glycan_intensity_dic):
        """
        Generate a set of linearized glycan sequences, from glycan_intensity_dictionary.
        Here we use a dictionary to store the results of de novo sequencing.
        Time complexity should be O(n^4)
        :param glycan: the number representation of glycans, such as "5410"
        :param glycan_intensity_dic: a dictionary contains composition and intensity
            Such as {'2410': 400.3, '1300': 589.3, '2310': 491.9, '4300': 413.3, '0100': 3591.5,
            '4310': 465.1, '1200': 648.3, '2200': 1793.2, '3200': 1706.6}
        :return: original amino acid name string
        """

        # Check the input types
        if not isinstance(glycan, (str,)):
            raise ValueError("glycan must be a string")
        if not isinstance(glycan_intensity_dic, (dict,)):
            raise ValueError("glycan_precursor must be a dictionary")

        # First step, count total number of glycan, and generate a list with the length of it.
        glycan_length = self.count_total(glycan)
        # total monosaccharides, get the maximum ions sequence number, index [0] will be an empty list here.
        glycan_lists = [[] for _ in range(glycan_length)]

        # Second step, count the total number of monosaccharides for the ions
        for key, value in glycan_intensity_dic.items():
            ion_number = self.count_total(key)
            ion_dic = {key: value}
            glycan_lists[ion_number].append(ion_dic)

        # Third step, generate glycan linear set based on the former ones
        glycan_dic = {}

        for i in range(1, glycan_length):
            for j in range(len(glycan_lists[i])):
                # Base case for only one monosaccharide
                if i == 1:
                    for key, value in glycan_lists[i][j].items():
                        monosaccharide = monosaccharide_number_replacement.get(key)
                        glycan_dic[monosaccharide] = value
                    print(glycan_dic)
                # Dynamic programming for generating ladders by top-down methods, from Y-ions to B-ions
                else:
                    # Generate all the possible candidates, based on former results in the set.
                    for key, value in glycan_lists[i][j].items():
                        for former_key in list(glycan_dic):
                            former_glycan = self.glycan_list_to_str(former_key)
                            if self.glycan_contain(former_glycan, key):
                                glycan_remains = self.glycan_minus(former_glycan, key)
                                print(glycan_remains)
                                glycan_dic_lists = self.glycan_str_to_list(glycan_remains)
                                for glycan_dic_list in glycan_dic_lists:
                                    print(former_key)
                                    new_key = former_key + "".join(glycan_dic_list)
                                    print(new_key)
                                    print(value)
                                    new_value = round(glycan_dic[former_key] + value, 1)
                                    glycan_dic[new_key] = new_value
                                    print(glycan_dic)

        for key in list(glycan_dic):
            print(key + ": " + str(glycan_dic[key]))

        return glycan_dic

