import time
import os

from PLNE_for_Delta1m_m1_decomposition import decompose as decompose_under_L3

from deformationAlgorithmForL3WithConstraintAboveCriteriaOrder import \
    recommendation_after_deformation_L3 as translationL3

TRANSLATION_METHODS = [translationL3]
from datetime import datetime

m = 9  # a lancer
A_dest_directory = f'./DONNEES/m_{m}/A_directory_m{m}'
W_dest_directory = f'./DONNEES/m_{m}/W_directory_m{m}'


def recommendation_algorithm(AlternativesSubsetsList, Wdict, decomposition_function=None):
    SortedAlternativesSubsetsList = sorted(AlternativesSubsetsList,
                                           key=lambda alt: sum([Wdict[criterion] for criterion in alt]), reverse=True)
    S_ = [None]
    for v in range(1, len(SortedAlternativesSubsetsList)):
        proSet, conSet = SortedAlternativesSubsetsList[0], \
                         SortedAlternativesSubsetsList[v]
        success_v, Sv = decomposition_function(proSet, conSet, Wdict)
        if not success_v:
            return False, None
        else:
            S_.append(Sv)

    return True, S_


def recommendation_algorithm_any_tree_height(AlternativesSubsetsList, W_dict, decomposition_function=None):
    SortedAlternativesSubsetsList = sorted(AlternativesSubsetsList,
                                           key=lambda alt: sum([W_dict[criterion] for criterion in alt]), reverse=True)
    SS = dict()
    for v in range(1, len(SortedAlternativesSubsetsList)):
        success_v = False
        for u in range(0, v):
            proSet, conSet = SortedAlternativesSubsetsList[u], \
                             SortedAlternativesSubsetsList[v]
            success_uv, Suv = decomposition_function(proSet, conSet, W_dict)
            if success_uv:
                success_v = True
                SS[(u, v)] = Suv
                break
        if not success_v:
            return False, None

    return True, SS


def recommendation_relaxation(AlternativesSubsetsList, Wdict, decomposition_function=None, kmax=2):
    SortedAlternativesSubsetsList = sorted(AlternativesSubsetsList,
                                           key=lambda alt: sum([Wdict[criterion] for criterion in alt]), reverse=True)
    # print(SortedAlternativesSubsetsList)

    S = dict()
    k = 2
    while k <= kmax + 1:
        v = k
        while v <= len(SortedAlternativesSubsetsList):
            S[v] = dict()
            u = 1
            failure_uv = False
            while not failure_uv and u <= k - 1:
                proSet, conSet = SortedAlternativesSubsetsList[u - 1], \
                                 SortedAlternativesSubsetsList[v - 1]
                success_uv, Suv = decomposition_function(proSet, conSet, Wdict)
                if not success_uv:
                    failure_uv = True
                    break
                S[v][u] = Suv
                u += 1
            if failure_uv:
                k += 1
                break
            v += 1
        if v <= len(SortedAlternativesSubsetsList):
            del S[v]
            continue

        return True, S

    return False, None


def check_int_List_1(int_List):
    A_binary_encoding_list = [integer_to_bin(val_int) for val_int in int_List]
    """significativite des m criteres"""
    for i in range(m):
        if all([alt[i] == '0' for alt in A_binary_encoding_list]) or all(
                [alt[i] == '1' for alt in A_binary_encoding_list]):
            return False

    return True


def check_int_List_2(int_List):
    A_binary_encoding_list = [integer_to_bin(val_int) for val_int in int_List]
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


def integer_to_bin(val_int):
    return format(val_int, "b").zfill(m)


def load_set_of_alternatives_set(filename):
    Alternatives_set_list = list()
    with open(filename) as alt_file:
        n = int(alt_file.readline())
        line = alt_file.readline()
        while line != "":
            line = line[:len(line) - 1]  # supprimer '\n' de fin
            line_alt_int_List = [integer_to_bin(int(e)) for e in line.split(" ", maxsplit=n - 1)]
            line_alt_subset_List = [{chr(ord('a') + i) for i in range(m) if alt[i] == '1'} for alt in line_alt_int_List]
            Alternatives_set_list.append(line_alt_subset_List)
            line = alt_file.readline()

    return Alternatives_set_list


def load_set_of_score_functions(filename):
    Omega_list = list()
    with open(filename) as w_file:
        line = w_file.readline()
        while line != "":
            line = line[:len(line) - 1]  # supprimer '\n' de fin
            # print(line.split(' ', maxsplit=m-1))
            line_score_function = [float(e) for e in line.split(' ', maxsplit=m - 1)]
            score_function_dict = {chr(ord('a') + i): line_score_function[i] for i in range(m)}
            Omega_list.append(score_function_dict)
            line = w_file.readline()
        return Omega_list


if __name__ == "__main__":
    LANGUAGES = [decompose_under_L3]
    assert len(LANGUAGES) == 1

    for A_list_file in os.listdir(A_dest_directory):
        for W_list_file in os.listdir(W_dest_directory):
            A_list = load_set_of_alternatives_set(A_dest_directory + '/' + A_list_file)
            W_list = load_set_of_score_functions(W_dest_directory + '/' + W_list_file)
            print("\n", A_list_file, "||", W_list_file)

            RESULT_RECO_ETOILE = {i: 0 for i in
                                  range(0, len(LANGUAGES))}
            RESULT_RECO_TRANSITIVITY = {i: 0 for i in range(0, len(LANGUAGES))}
            RESULT_RECO_RELAXATION = {i: 0 for i in range(0, len(LANGUAGES))}
            RESULT_RECO_TRANSLATION = {i: 0 for i in range(0, len(LANGUAGES))}
            print("Date", datetime.now())
            number = 0
            tmps_1, tmps_N, tmps_relax, tmps_translation = 0, 0, 0, 0

            for indice in range(len(A_list)):
                A = A_list[indice]
                W = W_list[indice]

                if number % 500 == 0:
                    print(number, datetime.now())
                    if number != 0:
                        print("*:", RESULT_RECO_ETOILE[0] / number, "tree:", RESULT_RECO_TRANSITIVITY[0] / number,
                              "relax:", RESULT_RECO_RELAXATION[0] / number, "translation:",
                              RESULT_RECO_TRANSLATION[0] / number)
                number += 1
                CPT_DICT = {i: False for i in range(0, len(LANGUAGES))}
                ilanguage = 0
                language = LANGUAGES[ilanguage]
                success_h1, S_h1 = recommendation_algorithm(A, W, decomposition_function=language)
                if success_h1:
                    RESULT_RECO_ETOILE[0] += 1
                    RESULT_RECO_TRANSITIVITY[0] += 1
                    RESULT_RECO_RELAXATION[0] += 1
                    RESULT_RECO_TRANSLATION[0] += 1
                else:
                    deb = time.time()
                    success_hN, S_hN = recommendation_algorithm_any_tree_height(A, W,
                                                                                decomposition_function=language)
                    if success_hN:
                        RESULT_RECO_TRANSITIVITY[0] += 1

                    success_relax, S_relax = recommendation_relaxation(A, W, decomposition_function=language)
                    if success_relax:
                        RESULT_RECO_RELAXATION[0] += 1

                    success_translation, S_translation = TRANSLATION_METHODS[ilanguage](A, W)
                    if success_translation:
                        RESULT_RECO_TRANSLATION[0] += 1

            print("RECO-ETOILE",
                  [100. * RESULT_RECO_ETOILE[il] / len(A_list)  for il in range(0, len(LANGUAGES))],
                  "%", "\tDurée :", tmps_1 / len(A_list) , "secondes")

            print("RECO-TRANSITIVITY", [100. * RESULT_RECO_TRANSITIVITY[il] / len(A_list) for il in
                                        range(0, len(LANGUAGES))], "%", "\tDurée :", tmps_N / len(A_list),
                  "secondes")

            print("RECO-RELAXATION", [100. * RESULT_RECO_RELAXATION[il] / len(A_list) for il in
                                      range(0, len(LANGUAGES))], "%", "\tDurée :",
                  tmps_relax / len(A_list), "secondes")

            print("RECO-TRANSLATION", [100. * RESULT_RECO_TRANSLATION[il] / len(A_list) for il in
                                       range(0, len(LANGUAGES))], "%", "\tDurée :",
                  tmps_translation / len(A_list), "secondes")
    print("Date", datetime.now())
