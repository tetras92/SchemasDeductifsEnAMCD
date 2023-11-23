import numpy as np


def integer_to_bin(val_int, m):
    return format(val_int, "b").zfill(m)


def check_int_List_1(int_List, m):
    A_binary_encoding_list = [integer_to_bin(val_int, m) for val_int in int_List]
    """significativite des m criteres"""
    for i in range(m):
        if all([alt[i] == '0' for alt in A_binary_encoding_list]) or all(
                [alt[i] == '1' for alt in A_binary_encoding_list]):
            return False

    return True


def check_int_List_2(int_List, m):
    A_binary_encoding_list = [integer_to_bin(val_int, m) for val_int in int_List]
    """absence de dominance"""
    for alt1 in A_binary_encoding_list:
        for alt2 in A_binary_encoding_list:
            if alt2 != alt1:
                alt1_int = [int(alt1[i]) for i in range(m)]
                alt2_int = [int(alt2[i]) for i in range(m)]
                if all([alt1_int[i] >= alt2_int[i] for i in range(m)]) or all(
                        [alt2_int[i] >= alt1_int[i] for i in range(m)]):
                    return False
    return True


# print([int(dig) for dig in integer_to_bin(4, 4)])

def predicat_in_for_arrays(arr, List_of_arrays):
    for oth_arr in List_of_arrays:
        if all(np.array(arr) == np.array(oth_arr)):
            return True
    return False


def save_instance_ijar(repository, filename, omega_score, setA_as_binary_np_arrays, setAR_as_binary_np_arrays,
                       AddedInfo):
    PI_set, Oth_set, Subj_set = AddedInfo
    with open(repository + filename, 'w') as ijarfile:
        ijarfile.write("w\n")
        str_omega_score = ""
        for val in omega_score:
            str_omega_score += str(val) + " "
        ijarfile.write(str_omega_score + "\n")
        ijarfile.write("A\n")
        for k in range(len(setA_as_binary_np_arrays)):
            elmt = setA_as_binary_np_arrays[k]
            str_elmt = str(elmt)
            if k in PI_set:
                str_elmt += " p"
            elif k in Oth_set and k not in Subj_set:
                str_elmt += " o"
            elif k in Oth_set and k in Subj_set:
                str_elmt += " os"
            str_elmt += "\n"
            ijarfile.write(str_elmt)

        ijarfile.write("AR\n")
        for elmt in setAR_as_binary_np_arrays:
            ijarfile.write(str(elmt) + "\n")
    ijarfile.close()


def save_instance_ijar2(repository, filename, omega_score, setA_as_binary_np_arrays, setAR_as_binary_np_arrays,
                        AddedInfo):
    PI_set, Oth_set, Subj_set, subj_othSet = AddedInfo
    with open(repository + filename, 'w') as ijarfile:
        ijarfile.write("w\n")
        str_omega_score = ""
        for val in omega_score:
            str_omega_score += str(val) + " "
        ijarfile.write(str_omega_score + "\n")
        ijarfile.write("A\n")
        for k in range(len(setA_as_binary_np_arrays)):
            elmt = setA_as_binary_np_arrays[k]
            str_elmt = str(elmt)
            if k in PI_set:
                str_elmt += "\tp"
            elif k in Oth_set and k not in Subj_set:
                str_elmt += "\to"  # remettre o pour indiquer que non eligible
            elif k in Oth_set and k not in subj_othSet:
                str_elmt += "\tx"
            elif k in Oth_set and k in subj_othSet:
                str_elmt += "\t"
                for pair in subj_othSet[k]:
                    str_elmt += str(pair) + " "
            str_elmt += "\n"
            ijarfile.write(str_elmt)

        ijarfile.write("AR\n")
        for elmt in setAR_as_binary_np_arrays:
            ijarfile.write(str(elmt) + "\n")
    ijarfile.close()


def save_instance_ijar3(repository, filename, omega_score, setA_as_binary_np_arrays, setAR_as_binary_np_arrays,
                        AddedInfo):
    PI_set, Oth_set, Subj_set, subj_othSet = AddedInfo
    with open(repository + filename, 'w') as ijarfile:
        ijarfile.write("w\n")
        str_omega_score = ""
        for val in omega_score:
            str_omega_score += str(val) + " "
        ijarfile.write(str_omega_score + "\n")
        ijarfile.write("A\n")
        for k in range(len(setA_as_binary_np_arrays)):
            elmt = setA_as_binary_np_arrays[k]
            str_elmt = str(elmt)
            if k in PI_set:
                str_elmt += "\tp"
            elif k in Oth_set and k not in Subj_set:
                str_elmt += "\to"  # remettre o pour indiquer que non eligible
            elif k in Oth_set and k not in subj_othSet:
                str_elmt += "\tx"
            elif k in Oth_set and k in subj_othSet:
                str_elmt += "\t"
                for pair in subj_othSet[k][0]:  # necessary
                    # str_elmt += f'*{str(pair)} '
                    str_elmt += f'*{pair[0]}>{pair[1]} '
                for pair in subj_othSet[k][1]:  # not necessary
                    # str_elmt += f'~{str(pair)} '
                    str_elmt += f'~{pair[0]}>{pair[1]} '
            str_elmt += "\n"
            ijarfile.write(str_elmt)

        ijarfile.write("AR\n")
        for elmt in setAR_as_binary_np_arrays:
            ijarfile.write(str(elmt) + "\n")
    ijarfile.close()


def save_instance_ijar4(repository, filename, omega_score, setA_as_binary_np_arrays, setAR_as_binary_np_arrays,
                        AddedInfo):
    PI_set, Oth_set, Subj_set, subj_othSet = AddedInfo
    with open(repository + filename, 'w') as ijarfile:
        ijarfile.write("w\n")
        str_omega_score = ""
        for val in omega_score:
            str_omega_score += str(val) + " "
        ijarfile.write(str_omega_score + "\n")
        ijarfile.write("A\n")
        for k in range(len(setA_as_binary_np_arrays)):
            elmt = setA_as_binary_np_arrays[k]
            str_elmt = str(elmt)
            if k in PI_set:
                str_elmt += "\tp"
            elif k in Oth_set and k not in Subj_set:
                str_elmt += "\to"  # remettre o pour indiquer que non eligible
            elif k in Oth_set and k not in subj_othSet:
                str_elmt += "\tx"
            elif k in Oth_set and k in subj_othSet:
                str_elmt += "\t"
                if len(subj_othSet[k]) == 1:
                    for pair in subj_othSet[k]["no-side"][0]:  # necessary
                        str_elmt += f'*{pair[0]}>{pair[1]} '
                    for pair in subj_othSet[k]["no-side"][1]:  # not necessary
                        str_elmt += f'~{pair[0]}>{pair[1]} '
                else:
                    for pair in subj_othSet[k]["side-x"][0]:  # necessary
                        str_elmt += f'*{pair[0]}>{pair[1]} '
                    for pair in subj_othSet[k]["side-x"][1]:  # not necessary
                        str_elmt += f'~{pair[0]}>{pair[1]} '
                    # str_elmt += f'# ({str(subj_othSet[k]["selected"][0])}) {str(subj_othSet[k]["selected"][1])}->{str(subj_othSet[k]["selected"][3])} *{str(subj_othSet[k]["selected"][2])}* # '
                    str_elmt += f'#({str(subj_othSet[k]["selected"][0])})->({str(subj_othSet[k]["selected"][2])})# '
                    for pair in subj_othSet[k]["side-y"][0]:  # necessary
                        str_elmt += f'*{pair[0]}>{pair[1]} '
                    for pair in subj_othSet[k]["side-y"][1]:  # not necessary
                        str_elmt += f'~{pair[0]}>{pair[1]} '
            str_elmt += "\n"
            ijarfile.write(str_elmt)

        ijarfile.write("AR\n")
        for elmt in setAR_as_binary_np_arrays:
            ijarfile.write(str(elmt) + "\n")
    ijarfile.close()


def save_instance_ijar6(repository, filename, omega_score, setA_as_binary_np_arrays, setAR_as_binary_np_arrays,
                        AddedInfo):
    PI_set, Oth_set, Subj_set, subj_othSet = AddedInfo
    with open(repository + filename, 'w') as ijarfile:
        ijarfile.write("w\n")
        str_omega_score = ""
        for val in omega_score:
            str_omega_score += str(val) + " "
        ijarfile.write(str_omega_score + "\n")
        ijarfile.write("A\n")

        for k in range(len(setA_as_binary_np_arrays) // 2):
            elmt = setA_as_binary_np_arrays[k]
            str_elmt = str(elmt)
            str_elmt += '\n'
            ijarfile.write(str_elmt)
        for k in range(len(setA_as_binary_np_arrays) // 2, len(setA_as_binary_np_arrays)):
            elmt = setA_as_binary_np_arrays[k]
            str_elmt = str(elmt) + "\t@"
            for i in range(len(setA_as_binary_np_arrays) // 2):
                if (i, k) in PI_set:
                    str_elmt += "\tp\t"
                elif (i, k) in Oth_set and (i, k) not in Subj_set:
                    str_elmt += "\to\t"  # remettre o pour indiquer que non eligible
                elif (i, k) in Oth_set and (i, k) not in subj_othSet:
                    str_elmt += "\tx"
                elif (i, k) in Oth_set and (i, k) in subj_othSet:
                    Liste_swaps_necessaires, Liste_swaps_non_necessaires = subj_othSet[(i, k)]
                    for pair in Liste_swaps_necessaires:
                        str_elmt += f'*{pair[0]}>{pair[1]} '
                    for pair in Liste_swaps_non_necessaires:
                        str_elmt += f'~{pair[0]}>{pair[1]} '
                else:
                    str_elmt += "\t-"
                str_elmt += " \t@"
            str_elmt += "\n"
            ijarfile.write(str_elmt)

        ijarfile.write("AR\n")
        for elmt in setAR_as_binary_np_arrays:
            ijarfile.write(str(elmt) + "\n")
    ijarfile.close()


def save_instance_ijar6a(repository, filename, omega_score, setA_as_binary_np_arrays, setAR_as_binary_np_arrays,
                         AddedInfo):
    PI_set, Oth_set, Subj_set, subj_othSet, dict_support = AddedInfo
    with open(repository + filename, 'w') as ijarfile:
        ijarfile.write("w\n")
        str_omega_score = ""
        for val in omega_score:
            str_omega_score += str(val) + " "
        ijarfile.write(str_omega_score + "\n")
        ijarfile.write("A\n")

        for k in range(len(setA_as_binary_np_arrays) // 2):
            elmt = setA_as_binary_np_arrays[k]
            str_elmt = str(elmt)
            str_elmt += '\n'
            ijarfile.write(str_elmt)
        for k in range(len(setA_as_binary_np_arrays) // 2, len(setA_as_binary_np_arrays)):
            elmt = setA_as_binary_np_arrays[k]
            str_elmt = str(elmt) + "\t@"
            for i in range(len(setA_as_binary_np_arrays) // 2):
                if (i, k) in PI_set:
                    str_elmt += "\tp\t"
                elif (i, k) in Oth_set and (i, k) not in Subj_set:
                    str_elmt += "\to\t"  # remettre o pour indiquer que non eligible
                elif (i, k) in Oth_set and (i, k) not in subj_othSet:
                    str_elmt += "\tx"
                elif (i, k) in Oth_set and (i, k) in subj_othSet:
                    Liste_swaps_necessaires, Liste_swaps_non_necessaires = subj_othSet[(i, k)]
                    for pair in Liste_swaps_necessaires:
                        str_elmt += f'*{pair[0]}>{pair[1]} '
                    for pair in Liste_swaps_non_necessaires:
                        str_elmt += f'~{pair[0]}>{pair[1]} '
                    str_elmt += f'. {{{dict_support[(i, k)]} %}} '
                else:
                    str_elmt += "\t-"
                str_elmt += " \t@"
            str_elmt += "\n"
            ijarfile.write(str_elmt)

        ijarfile.write("AR\n")
        for elmt in setAR_as_binary_np_arrays:
            ijarfile.write(str(elmt) + "\n")
    ijarfile.close()


def read_ijar6(repository, filename):
    omega_score_function = None
    PI_set, Oth_set, Subj_set, subj_othSet_dict = set(), set(), set(), dict()
    setA_as_binary_np_arrays = list()
    setAR_as_binary_np_arrays = list()
    with open(repository + filename, 'r') as ijarfile:
        _ = ijarfile.readline()
        line_omega_score = ijarfile.readline()
        omega_score_function = np.array([float(x) for x in line_omega_score.split(' ')[:-1]])
        _ = ijarfile.readline()  # a
        k_id_line = 0
        while True:
            line = ijarfile.readline()
            if line == "AR\n":
                break
            arr, addinfo = line.split(']')
            arr = arr[1:]  # retirer [
            setA_as_binary_np_arrays.append(np.array([int(x) for x in arr.split(' ')]))
            if not (addinfo == "\n"):
                addinfo = addinfo[2:-2]  # retirer l'espace, premier @, dernier @ et \n
                addinfoList = addinfo.split("@")
                for i_ in range(len(addinfoList)):
                    xpl_i_k = addinfoList[i_].strip()  # retirer espace et tabulation autour
                    if xpl_i_k == "p":
                        PI_set.add((i_, k_id_line))
                    elif xpl_i_k == "o" or xpl_i_k == "x":
                        Oth_set.add((i_, k_id_line))
                    elif xpl_i_k == "-":
                        pass
                    else:
                        Oth_set.add((i_, k_id_line))
                        Subj_set.add((i_, k_id_line))
                        explanation_str_list = xpl_i_k.split(' ')
                        necessary_swaps_list, non_necessary_swaps_list = list(), list()
                        for chainons_str in explanation_str_list:
                            chainons_str_without_prefix = chainons_str[1:]
                            chainon_pair = chainons_str_without_prefix.split('>')
                            if chainons_str[0] == "*":
                                necessary_swaps_list.append((int(chainon_pair[0]), int(chainon_pair[1])))
                            elif chainons_str[0] == "~":
                                non_necessary_swaps_list.append((int(chainon_pair[0]), int(chainon_pair[1])))
                            else:
                                raise Exception("prefix epistemiq non connue : ", "a", chainons_str[0])
                        subj_othSet_dict[(i_, k_id_line)] = (necessary_swaps_list, non_necessary_swaps_list)
            k_id_line += 1

        line = ijarfile.readline()
        while line != "":
            line = line[1:-2]  # retirer [ et ]\n
            setAR_as_binary_np_arrays.append(np.array([int(x) for x in line.split(' ')]))
            line = ijarfile.readline()
    ijarfile.close()

    return omega_score_function, setA_as_binary_np_arrays, setAR_as_binary_np_arrays, (
        PI_set, Oth_set, Subj_set, subj_othSet_dict)


def read_ijar6a(repository, filename):
    omega_score_function = None
    PI_set, Oth_set, Subj_set, subj_othSet_dict, support_Dict = set(), set(), set(), dict(), dict()
    setA_as_binary_np_arrays = list()
    setAR_as_binary_np_arrays = list()
    with open(repository + filename, 'r') as ijarfile:
        _ = ijarfile.readline()
        line_omega_score = ijarfile.readline()
        omega_score_function = np.array([float(x) for x in line_omega_score.split(' ')[:-1]])
        _ = ijarfile.readline()  # a
        k_id_line = 0
        while True:
            line = ijarfile.readline()
            if line == "AR\n":
                break
            arr, addinfo = line.split(']')
            arr = arr[1:]  # retirer [
            setA_as_binary_np_arrays.append(np.array([int(x) for x in arr.split(' ')]))
            if not (addinfo == "\n"):
                addinfo = addinfo[2:-2]  # retirer l'espace, premier @, dernier @ et \n
                addinfoList = addinfo.split("@")
                for i_ in range(len(addinfoList)):
                    xpl_i_k_and_support = addinfoList[i_].strip()  # retirer espace et tabulation autour
                    if xpl_i_k_and_support == "p":
                        PI_set.add((i_, k_id_line))
                    elif xpl_i_k_and_support == "o" or xpl_i_k_and_support == "x":
                        Oth_set.add((i_, k_id_line))
                    elif xpl_i_k_and_support == "-":
                        pass
                    else:
                        Oth_set.add((i_, k_id_line))
                        Subj_set.add((i_, k_id_line))
                        part_xp, part_support = xpl_i_k_and_support.split(".")
                        part_xp = part_xp.strip()
                        part_support = part_support.strip()
                        explanation_str_list = part_xp.split(' ')
                        necessary_swaps_list, non_necessary_swaps_list = list(), list()
                        for chainons_str in explanation_str_list:
                            chainons_str_without_prefix = chainons_str[1:]
                            chainon_pair = chainons_str_without_prefix.split('>')
                            if chainons_str[0] == "*":
                                necessary_swaps_list.append((int(chainon_pair[0]), int(chainon_pair[1])))
                            elif chainons_str[0] == "~":
                                non_necessary_swaps_list.append((int(chainon_pair[0]), int(chainon_pair[1])))
                            else:
                                raise Exception("prefix epistemiq non connue : ", chainons_str[0])
                        subj_othSet_dict[(i_, k_id_line)] = (necessary_swaps_list, non_necessary_swaps_list)
                        support_Dict[(i_, k_id_line)] = int(part_support[1:].split(" ")[0])                     # bricolage
            k_id_line += 1

        line = ijarfile.readline()
        while line != "":
            line = line[1:-2]  # retirer [ et ]\n
            setAR_as_binary_np_arrays.append(np.array([int(x) for x in line.split(' ')]))
            line = ijarfile.readline()
    ijarfile.close()

    return omega_score_function, setA_as_binary_np_arrays, setAR_as_binary_np_arrays, (
        PI_set, Oth_set, Subj_set, subj_othSet_dict, support_Dict)


def save_instance_ijar5(repository, filename, omega_score, setA_as_binary_np_arrays, setAR_as_binary_np_arrays,
                        AddedInfo):
    new_PI_set, new_Oth_set, new_Deduced, new_Deduced_trivial, new_Deduced_non_trivial, new_Deduced_Subject, Details_Dict = AddedInfo
    with open(repository + filename, 'w') as ijarfile:
        ijarfile.write("w\n")
        str_omega_score = ""
        for val in omega_score:
            str_omega_score += str(val) + " "
        ijarfile.write(str_omega_score + "\n")
        ijarfile.write("A\n")

        for k in range(len(setA_as_binary_np_arrays) // 2):
            elmt = setA_as_binary_np_arrays[k]
            str_elmt = str(elmt)
            str_elmt += '\n'
            ijarfile.write(str_elmt)
        for k in range(len(setA_as_binary_np_arrays) // 2, len(setA_as_binary_np_arrays)):
            elmt = setA_as_binary_np_arrays[k]
            str_elmt = str(elmt) + "\t@ \t\t"
            for i in range(len(setA_as_binary_np_arrays) // 2):
                if (i, k) in new_PI_set:
                    str_elmt += "\tp\t"
                elif (i, k) in new_Deduced_trivial:
                    str_elmt += "\tt\t"
                elif (i, k) in new_Deduced_non_trivial and (i, k) not in new_Deduced_Subject:
                    str_elmt += "\to\t"  # remettre o pour indiquer que non eligible
                elif (i, k) in new_Deduced_Subject and (i, k) not in Details_Dict:
                    str_elmt += "\tx\t"
                elif (i, k) in Details_Dict:
                    str_elmt += "\t"
                    Liste_swaps_necessaires, Liste_swaps_non_necessaires, indexes = Details_Dict[(i, k)]
                    if Liste_swaps_non_necessaires is None:
                        if indexes is None:  # explication necessaire
                            for pair in Liste_swaps_necessaires:
                                str_elmt += f'*{pair[0]}>{pair[1]} '
                        else:  # generalized swap necessaires et PI
                            # index_alt_z1, index_alt_z2 = indexes                # Cas 1_generalized
                            # for pair in Liste_swaps_necessaires[0]:
                            #     str_elmt += f'*{pair[0]}>{pair[1]} '
                            # str_elmt += f'#({index_alt_z1})->({index_alt_z2})# '
                            # for pair in Liste_swaps_necessaires[1]:
                            #     str_elmt += f'*{pair[0]}>{pair[1]} '
                            for i_ in range(len(indexes)):
                                for pair in Liste_swaps_necessaires[i_]:
                                    str_elmt += f'*{pair[0]}>{pair[1]} '
                                index_alt_z1, index_alt_z2 = indexes[i_]
                                str_elmt += f'#({index_alt_z1})->({index_alt_z2})# '
                            for pair in Liste_swaps_necessaires[-1]:
                                str_elmt += f'*{pair[0]}>{pair[1]} '
                    elif indexes is None:  # mixed simple
                        for pair in Liste_swaps_necessaires:
                            str_elmt += f'*{pair[0]}>{pair[1]} '
                        for pair in Liste_swaps_non_necessaires:
                            str_elmt += f'~{pair[0]}>{pair[1]} '
                    else:  # mixed 1-generalized
                        index_alt_z1, index_alt_z2 = indexes
                        for pair in Liste_swaps_necessaires[0]:
                            str_elmt += f'*{pair[0]}>{pair[1]} '
                        for pair in Liste_swaps_non_necessaires[0]:
                            str_elmt += f'~{pair[0]}>{pair[1]} '
                        str_elmt += f'#({index_alt_z1})->({index_alt_z2})# '
                        for pair in Liste_swaps_necessaires[1]:
                            str_elmt += f'*{pair[0]}>{pair[1]} '
                        for pair in Liste_swaps_non_necessaires[1]:
                            str_elmt += f'~{pair[0]}>{pair[1]} '
                str_elmt += " \t@"
            str_elmt += "\n"
            ijarfile.write(str_elmt)

        ijarfile.write("AR\n")
        for elmt in setAR_as_binary_np_arrays:
            ijarfile.write(str(elmt) + "\n")
    ijarfile.close()


def read_ijar(repository, filename):
    omega_score_function = None
    PI_set, Oth_set, Subj_set = set(), set(), set()
    setA_as_binary_np_arrays = list()
    setAR_as_binary_np_arrays = list()
    with open(repository + filename, 'r') as ijarfile:
        str_w = ijarfile.readline()
        line_omega_score = ijarfile.readline()
        omega_score_function = np.array([float(x) for x in line_omega_score.split(' ')[:-1]])
        _ = ijarfile.readline()  # a
        id_line = 0
        while True:
            line = ijarfile.readline()
            if line == "AR\n":
                break
            arr, addinfo = line.split(']')
            arr = arr[1:]  # retirer [
            setA_as_binary_np_arrays.append(np.array([int(x) for x in arr.split(' ')]))
            if not (addinfo == "\n"):
                addinfo = addinfo[1:-1]  # retirer l'espace et \n
                if addinfo == "p":
                    PI_set.add(id_line)
                elif addinfo == "o":
                    Oth_set.add(id_line)
                elif addinfo == "os":
                    Oth_set.add(id_line)
                    Subj_set.add(id_line)
                else:
                    raise Exception("Chaine differente de p, o, os : " + addinfo)
            id_line += 1

        line = ijarfile.readline()
        while line != "":
            line = line[1:-2]  # retirer [ et ]\n
            setAR_as_binary_np_arrays.append(np.array([int(x) for x in line.split(' ')]))
            line = ijarfile.readline()
    ijarfile.close()

    return omega_score_function, setA_as_binary_np_arrays, setAR_as_binary_np_arrays, (PI_set, Oth_set, Subj_set)


def read_ijar3(repository, filename):
    omega_score_function = None
    PI_set, Oth_set, Subj_set, subj_othSet_dict = set(), set(), set(), dict()
    setA_as_binary_np_arrays = list()
    setAR_as_binary_np_arrays = list()
    with open(repository + filename, 'r') as ijarfile:
        str_w = ijarfile.readline()
        line_omega_score = ijarfile.readline()
        omega_score_function = np.array([float(x) for x in line_omega_score.split(' ')[:-1]])
        _ = ijarfile.readline()  # a
        id_line = 0
        while True:
            line = ijarfile.readline()
            if line == "AR\n":
                break
            arr, addinfo = line.split(']')
            arr = arr[1:]  # retirer [
            setA_as_binary_np_arrays.append(np.array([int(x) for x in arr.split(' ')]))
            if not (addinfo == "\n"):
                addinfo = addinfo[1:-1]  # retirer l'espace et \n
                if addinfo == "p":
                    PI_set.add(id_line)
                elif addinfo == "o" or addinfo == "x":
                    Oth_set.add(id_line)
                elif addinfo == "os":
                    raise Exception("Pas de 'os' dans un .ijar3")
                else:  # il manque clairement le remplissage de Subj_Set et Oth_Set : 07/01/2023
                    explanation_str_list = addinfo.split('\t')[0].split(' ')[:-1]
                    necessary_swaps_list, non_necessary_swaps_list = list(), list()
                    for chainons_str in explanation_str_list:
                        chainons_str_without_prefix = chainons_str[1:]
                        chainon_pair = chainons_str_without_prefix.split('>')
                        if chainons_str[0] == "*":
                            necessary_swaps_list.append((int(chainon_pair[0]), int(chainon_pair[1])))
                        elif chainons_str[0] == "~":
                            non_necessary_swaps_list.append((int(chainon_pair[0]), int(chainon_pair[1])))
                        else:
                            raise Exception("prefix epistemiq non connue : ", chainons_str[0])
                    subj_othSet_dict[id_line] = (necessary_swaps_list, non_necessary_swaps_list)
                    # print(explanation_str_list)
            id_line += 1

        line = ijarfile.readline()
        while line != "":
            line = line[1:-2]  # retirer [ et ]\n
            setAR_as_binary_np_arrays.append(np.array([int(x) for x in line.split(' ')]))
            line = ijarfile.readline()
    ijarfile.close()

    return omega_score_function, setA_as_binary_np_arrays, setAR_as_binary_np_arrays, (
        PI_set, Oth_set, Subj_set, subj_othSet_dict)
