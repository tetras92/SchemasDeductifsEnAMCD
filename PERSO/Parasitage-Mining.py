import csv
from Reference_Functions_lv import *

import pandas as pd
import datetime
from NPComputationViaCofactor import np_computation

m_value = 4
directory_CBTO = f'/home/manuel239/PycharmProjects/SchemasDeductifsEnAMCD/CBTO{m_value}/'


def pair_of_sets_to_pair_of_array(set_pair, m):
    elmt0 = np.zeros(m)
    elmt1 = np.zeros(m)
    for i in set_pair[0]:
        elmt0[i] = 1
    for i in set_pair[1]:
        elmt1[i] = 1
    return elmt0, elmt1

def critical_irreductible(CriticalList):
    Critical_statements = [pair_of_sets_to_pair_of_array(crit_pair, m_value) for crit_pair in CriticalList]
    Critical_covectors = [st1 - st2 for st1, st2 in Critical_statements]
    Irreductibles = list()
    Reductibles = list()
    for k in range(len(Critical_statements)):
        OtherCovectors = Critical_covectors[0:k] + (
            Critical_covectors[k + 1:] if k < len(CriticalList) - 1 else [])
        OtherStatements = Critical_statements[0:k] + (
            Critical_statements[k + 1:] if k < len(CriticalList) - 1 else [])
        if not np_computation(Critical_covectors[k], OtherCovectors):
            Irreductibles.append(CriticalList[k])
        else:
            Reductibles.append(CriticalList[k])
    return Irreductibles

if __name__ == "__main__":
    ALL_SCHEMES_USED = set()
    DCompta = dict()
    DictAffine = dict()
    omega_score_function = None
    print(datetime.datetime.now())
    # L = list()
    min_cons = 1
    # max_cons = 100000
    max_cons = len(os.listdir(directory_CBTO))
    for id_model in range(1, len(os.listdir(directory_CBTO)) + 1):
        bool_tr = False
        bool_1_cg = False
        if id_model > max_cons or id_model < min_cons:
            continue
        if id_model != 4:
            continue
        omega_file = f'model{id_model}.csv'
        if id_model % 400 == 0:
            print(omega_file, datetime.datetime.now(), ALL_SCHEMES_USED)
        with open(directory_CBTO + omega_file) as utilityFile:
            reader = csv.DictReader(utilityFile)
            w_list = list()
            for row in reader:
                for criterion in reader.fieldnames:
                    w_list.append(float(row[criterion]))
        DCompta[omega_file] = dict()

        omega_score_function = np.array(w_list)

        omega_score_function_dict = {i: omega_score_function[i] for i in range(len(omega_score_function))}
        # if omega_score_function_dict != {0: 5, 1: 8, 2: 14, 3: 16, 4: 20}: continue
        # omega_score_function_dict =
        Critical_Pairs, AllDisjointsPairs = generate_critical_pairs_and_Tm(
            directory_CBTO + omega_file, m_value)
        Critical_Pairs = [({1,2}, {0,3}), ({0,2}, {1}), ({3}, {0,2})]
        print("=> model n°", id_model, ':', {chr(ord('a') + i): val for i, val in omega_score_function_dict.items()})
        print("=> Critical pairs", [("".join([chr(ord('a') + i) for i in sorted(st1)]), "".join([chr(ord('a') + i) for i in sorted(st2)])) for st1, st2 in Critical_Pairs])
        LesIrreductibles = critical_irreductible(Critical_Pairs)
        print("=> Indeductible Critical pairs", [("".join([chr(ord('a') + i) for i in sorted(st1)]), "".join([chr(ord('a') + i) for i in sorted(st2)])) for st1, st2 in LesIrreductibles])

        AllPairWiseComparisonsIndices = {k: ("?", None) for k in range(len(AllDisjointsPairs))}
        CorrespondancePWC_PSTR = {k: None for k in range(len(AllDisjointsPairs))}

        DictAffine[omega_file] = {"Base": np.NAN, "Compo2": np.NAN, "Compo3": np.NAN, "UNXP": np.NAN,
                                  "DynEvo": np.NAN,
                                  "NB_TO_EXPLAIN": len(AllDisjointsPairs), "NB_STEPS": np.NAN}
        # print(AllPairWiseComparisonsIndices)

        PI_statements = [pair_of_sets_to_pair_of_array(crit_pair, m_value) for crit_pair in Critical_Pairs]
        PI_cofactors = [st1 - st2 for st1, st2 in PI_statements]

        Xplained_statements = list()
        Xplained_cofactors = list()

        PairwiseIndicesXplained = set()
        step = 0
        stabilised = False
        Schemes_used = set()
        BDB = {"[tr]": dict(), "[cg]": dict(), "[cov]": dict(),
               "[cg/tr]": dict(), "[cg/cov]": dict(),
               "[tr/cg]": dict(), "[tr/cov]": dict(),
               "[cov/cg]": dict(), "[cov/tr]": dict(),
               "[cg/tr/cov]": dict(), "[cg/cov/tr]": dict(),
               "[tr/cg/cov]": dict(), "[tr/cov/cg]": dict(),
               "[cov/cg/tr]": dict(), "[cov/tr/cg]": dict(),
               }
        print("\n\n Let's explain \n")
        while not stabilised:
            stabilised = True
            TmpXplained_statements = list()
            TmpXplained_cofactors = list()
            DCompta[omega_file][step] = {"[tr]": 0, "[cg]": 0, "[cov]": 0,
                                         "[cg/tr]": 0, "[cg/cov]": 0,
                                         "[tr/cg]": 0, "[tr/cov]": 0,
                                         "[cov/cg]": 0, "[cov/tr]": 0,
                                         "[cg/tr/cov]": 0, "[cg/cov/tr]": 0,
                                         "[tr/cg/cov]": 0, "[tr/cov/cg]": 0,
                                         "[cov/cg/tr]": 0, "[cov/tr/cg]": 0,
                                         "UNXP": 0, "Considered": 0
                                         }

            for k in range(len(AllDisjointsPairs)):
                if k in PairwiseIndicesXplained:
                    continue
                DCompta[omega_file][step]["Considered"] += 1
                pair_to_explain = AllDisjointsPairs[k]
                statement = pair_of_sets_to_pair_of_array(pair_to_explain, m_value)
                cofactor = statement[0] - statement[1]

                statement_str = ("".join([chr(ord('a') + i) for i in range(len(statement[0])) if statement[0][i] == 1]),
                                 "".join([chr(ord('a') + i) for i in range(len(statement[1])) if statement[1][i] == 1]))
                CorrespondancePWC_PSTR[k] = statement_str

                # if statement_str != ('cd', 'e'):
                #     continue
                is_tr, d_tr = tr(statement, PI_statements + Xplained_statements, omega_score_function_dict, step=step,
                                 BDB=BDB)

                is_cg, d_cg = cg(statement, PI_statements + Xplained_statements, omega_score_function_dict, step=step,
                                 BDB=BDB)
                if is_cg:
                    AllPairWiseComparisonsIndices[k] = ("[cg]", d_cg)
                    stabilised = False
                    PairwiseIndicesXplained.add(k)
                    TmpXplained_statements.append(statement)
                    TmpXplained_cofactors.append(cofactor)
                    DCompta[omega_file][step]["[cg]"] += 1
                    Schemes_used.add("[cg]")
                    print(statement_str, d_cg)
                    continue

                if is_tr:
                    AllPairWiseComparisonsIndices[k] = ("[tr]", d_tr)
                    stabilised = False
                    PairwiseIndicesXplained.add(k)
                    TmpXplained_statements.append(statement)
                    TmpXplained_cofactors.append(cofactor)
                    DCompta[omega_file][step]["[tr]"] += 1
                    Schemes_used.add("[tr]")
                    print(statement_str, d_tr)
                    continue

                is_cov, d_cov = cov(statement, PI_statements + Xplained_statements, omega_score_function_dict,
                                    step=step,
                                    BDB=BDB)
                if is_cov:
                    AllPairWiseComparisonsIndices[k] = ("[cov]", d_cov)
                    stabilised = False
                    PairwiseIndicesXplained.add(k)
                    TmpXplained_statements.append(statement)
                    TmpXplained_cofactors.append(cofactor)
                    DCompta[omega_file][step]["[cov]"] += 1
                    Schemes_used.add("[cov]")
                    print(statement_str, d_cov)
                    continue

                is_trCg, d_trCg = tr_cg(statement, PI_statements + Xplained_statements, omega_score_function_dict,
                                        step=step,
                                        BDB=BDB)
                if is_trCg:
                    AllPairWiseComparisonsIndices[k] = ("[tr/cg]", d_trCg)
                    stabilised = False
                    PairwiseIndicesXplained.add(k)
                    TmpXplained_statements.append(statement)
                    TmpXplained_cofactors.append(cofactor)
                    DCompta[omega_file][step]["[tr/cg]"] += 1
                    Schemes_used.add("[tr/cg]")
                    continue

                is_trCov, d_trCov = tr_cov(statement, PI_statements + Xplained_statements, omega_score_function_dict,
                                           step=step,
                                           BDB=BDB)
                if is_trCov:
                    AllPairWiseComparisonsIndices[k] = ("[tr/cov]", d_trCov)
                    stabilised = False
                    PairwiseIndicesXplained.add(k)
                    TmpXplained_statements.append(statement)
                    TmpXplained_cofactors.append(cofactor)
                    DCompta[omega_file][step]["[tr/cov]"] += 1
                    Schemes_used.add("[tr/cov]")

                    continue
                is_cgTr, d_cgTr = cg_tr(statement, PI_statements + Xplained_statements, omega_score_function_dict,
                                        step=step,
                                        BDB=BDB)
                if is_cgTr:
                    AllPairWiseComparisonsIndices[k] = ("[cg/tr]", d_cgTr)
                    stabilised = False
                    PairwiseIndicesXplained.add(k)
                    TmpXplained_statements.append(statement)
                    TmpXplained_cofactors.append(cofactor)
                    DCompta[omega_file][step]["[cg/tr]"] += 1
                    Schemes_used.add("[cg/tr]")
                    print(statement_str, d_cgTr)
                    continue

                is_cgCov, d_cgCov = cg_cov(statement, PI_statements + Xplained_statements, omega_score_function_dict,
                                           step=step,
                                           BDB=BDB)
                if is_cgCov:
                    AllPairWiseComparisonsIndices[k] = ("[cg/cov]", d_cgCov)
                    stabilised = False
                    PairwiseIndicesXplained.add(k)
                    TmpXplained_statements.append(statement)
                    TmpXplained_cofactors.append(cofactor)
                    DCompta[omega_file][step]["[cg/cov]"] += 1
                    Schemes_used.add("[cg/cov]")
                    continue

                is_covCg, d_covCg = cov_cg(statement, PI_statements + Xplained_statements, omega_score_function_dict,
                                           step=step,
                                           BDB=BDB)
                if is_covCg:
                    AllPairWiseComparisonsIndices[k] = ("[cov/cg]", d_covCg)
                    stabilised = False
                    PairwiseIndicesXplained.add(k)
                    TmpXplained_statements.append(statement)
                    TmpXplained_cofactors.append(cofactor)
                    DCompta[omega_file][step]["[cov/cg]"] += 1
                    Schemes_used.add("[cov/cg]")
                    # print(statement_str, d_covCg)
                    continue

                is_covTr, d_covTr = cov_tr(statement, PI_statements + Xplained_statements, omega_score_function_dict,
                                           step=step,
                                           BDB=BDB)
                if is_covTr:
                    AllPairWiseComparisonsIndices[k] = ("[cov/tr]", d_covTr)
                    stabilised = False
                    PairwiseIndicesXplained.add(k)
                    TmpXplained_statements.append(statement)
                    TmpXplained_cofactors.append(cofactor)
                    DCompta[omega_file][step]["[cov/tr]"] += 1
                    Schemes_used.add("[cov/tr]")
                    # print(statement_str, d_covTr, id_model)
                    continue

                is_trCgCov, d_trCgCov = tr_generique(statement, PI_statements + Xplained_statements,
                                                     omega_score_function_dict, cg_cov,
                                                     step=step,
                                                     BDB=BDB)
                if is_trCgCov:
                    AllPairWiseComparisonsIndices[k] = ("[tr/cg/cov]", d_trCgCov)
                    stabilised = False
                    PairwiseIndicesXplained.add(k)
                    TmpXplained_statements.append(statement)
                    TmpXplained_cofactors.append(cofactor)
                    DCompta[omega_file][step]["[tr/cg/cov]"] += 1
                    Schemes_used.add("[tr/cg/cov]")
                    print(statement_str, d_trCgCov)
                    continue

                is_trCovCg, d_trCovCg = tr_generique(statement, PI_statements + Xplained_statements,
                                                     omega_score_function_dict, cov_cg,
                                                     step=step,
                                                     BDB=BDB)
                if is_trCovCg:
                    AllPairWiseComparisonsIndices[k] = ("[tr/cov/cg]", d_trCovCg)
                    stabilised = False
                    PairwiseIndicesXplained.add(k)
                    TmpXplained_statements.append(statement)
                    TmpXplained_cofactors.append(cofactor)
                    DCompta[omega_file][step]["[tr/cov/cg]"] += 1
                    Schemes_used.add("[tr/cov/cg]")
                    continue

                is_cgTrCov, d_cgTrCov = cg_generique(statement, PI_statements + Xplained_statements,
                                                     omega_score_function_dict, tr_cov,
                                                     step=step,
                                                     BDB=BDB)
                if is_cgTrCov:
                    AllPairWiseComparisonsIndices[k] = ("[cg/tr/cov]", d_cgTrCov)
                    stabilised = False
                    PairwiseIndicesXplained.add(k)
                    TmpXplained_statements.append(statement)
                    TmpXplained_cofactors.append(cofactor)
                    DCompta[omega_file][step]["[cg/tr/cov]"] += 1
                    Schemes_used.add("[cg/tr/cov]")
                    continue

                is_cgCovTr, d_cgCovTr = cg_generique(statement, PI_statements + Xplained_statements,
                                                     omega_score_function_dict, cov_tr,
                                                     step=step,
                                                     BDB=BDB)
                if is_cgCovTr:
                    AllPairWiseComparisonsIndices[k] = ("[cg/cov/tr]", d_cgCovTr)
                    stabilised = False
                    PairwiseIndicesXplained.add(k)
                    TmpXplained_statements.append(statement)
                    TmpXplained_cofactors.append(cofactor)
                    DCompta[omega_file][step]["[cg/cov/tr]"] += 1
                    Schemes_used.add("[cg/cov/tr]")
                    print(statement_str, d_cgCovTr)
                    continue

                is_covCgTr, d_covCgTr = cov_generique(statement, PI_statements + Xplained_statements,
                                                      omega_score_function_dict, cg_tr,
                                                      step=step,
                                                      BDB=BDB)
                if is_covCgTr:
                    AllPairWiseComparisonsIndices[k] = ("[cov/cg/tr]", d_covCgTr)
                    stabilised = False
                    PairwiseIndicesXplained.add(k)
                    TmpXplained_statements.append(statement)
                    TmpXplained_cofactors.append(cofactor)
                    DCompta[omega_file][step]["[cov/cg/tr]"] += 1
                    Schemes_used.add("[cov/cg/tr]")
                    continue

                is_covTrCg, d_covTrCg = cov_generique(statement, PI_statements + Xplained_statements,
                                                      omega_score_function_dict, tr_cg,
                                                      step=step,
                                                      BDB=BDB)
                if is_covTrCg:
                    AllPairWiseComparisonsIndices[k] = ("[cov/tr/cg]", d_covTrCg)
                    stabilised = False
                    PairwiseIndicesXplained.add(k)
                    TmpXplained_statements.append(statement)
                    TmpXplained_cofactors.append(cofactor)
                    DCompta[omega_file][step]["[cov/tr/cg]"] += 1
                    Schemes_used.add("[cov/tr/cg]")
                    # print(id_model)
                    # exit(239)
                    continue

                DCompta[omega_file][step]["UNXP"] += 1
                print(statement_str, "XXX")
            Xplained_statements += TmpXplained_statements
            Xplained_cofactors += TmpXplained_cofactors

            if not stabilised:
                # for k in AllPairWiseComparisonsIndices:
                #     print(CorrespondancePWC_PSTR[k], "<=====", AllPairWiseComparisonsIndices[k])
                # print(*AllPairWiseComparisonsIndices.items(), sep="\n")
                # print("Step =========>", step, f'({len(Schemes_used)})', Schemes_used)
                # print(datetime.datetime.now())
                To_remove = [(cat, elem) for cat, DCat in BDB.items() for elem, valueslist in DCat.items() if
                             not valueslist[0]]
                for cat, elem in To_remove:
                    del BDB[cat][elem]
            step += 1
            # break
        nb_total_xplained = 0
        ALL_SCHEMES_USED = ALL_SCHEMES_USED | Schemes_used
        L0, L1, L2, L3, C_L0, C_L1, C_L2, C_L3 = 0, 0, 0, 0, 0, 0, 0, 0
        AllStep0, OtherSteps = 0, 0
        UN_XP0 = None
        for stEp, det_stEp in DCompta[omega_file].items():
            BaseNumber = det_stEp['[tr]'] + det_stEp['[cg]'] + det_stEp['[cov]']
            Compo2Number = det_stEp['[cg/tr]'] + det_stEp['[cg/cov]'] + det_stEp["[tr/cg]"] + det_stEp["[tr/cov]"] + \
                           det_stEp["[cov/cg]"] + det_stEp["[cov/tr]"]
            Compo3Number = det_stEp['[cov/cg/tr]'] + det_stEp['[tr/cg/cov]'] + det_stEp["[cov/tr/cg]"] + det_stEp[
                "[cg/tr/cov]"] + det_stEp["[tr/cov/cg]"] + det_stEp["[cg/cov/tr]"]
            if stEp == 0:
                AllStep0 = BaseNumber + Compo2Number + Compo3Number
                DictAffine[omega_file]["Base"] = BaseNumber
                DictAffine[omega_file]["Compo2"] = Compo2Number
                DictAffine[omega_file]["Compo3"] = Compo3Number

                UN_XP = len(AllDisjointsPairs) - AllStep0
                DictAffine[omega_file]["UNXP"] = UN_XP
            else:
                OtherSteps += BaseNumber + Compo2Number + Compo3Number

            nb_total_xplained_stp = BaseNumber + Compo2Number + Compo3Number
            # print("Step", stEp, "Total XP-step : ", nb_total_xplained_stp)
            nb_total_xplained += nb_total_xplained_stp

        DictAffine[omega_file]["DynEvo"] = OtherSteps
        DictAffine[omega_file]["ALL_UNXP"] = len(AllDisjointsPairs) - nb_total_xplained
        DictAffine[omega_file]["NB_STEPS"] = step - 1

        if id_model % 300 == 0:
            print(datetime.datetime.now(), id_model)

    data = pd.DataFrame.from_dict(DictAffine, orient='index')
    # print(data.columns[0])
    # data.to_csv("STORAGE/mining-cov-" + str(min_cons) + "-" + str(max_cons) + ".csv")
    print(datetime.datetime.now())
    print(ALL_SCHEMES_USED)
    # print(L)
