import datetime

from NPComputationViaCovector import np_computation
import numpy as np
import csv
from itertools import combinations
from gurobipy import *
from Utils import save_instance_pfia, generate_instance
from X_robust_Engine import compute_robust_explanations

XP_DIR = f'./XP/'


def transform_omega_function_into_cbto(W_array):
    def ensemble_des_parties_sauf_vide(Ens):
        ALL = list()
        for i in range(1, len(Ens)):
            for elem in combinations(Ens, i):
                ALL.append(tuple(sorted(elem)))
        return ALL

    L = sorted(ensemble_des_parties_sauf_vide(set(range(len(W_array)))), key=lambda E: sum([W_array[e] for e in E]))
    # print(L[:len(L) // 2 + 1])

    modelcbto = Model("Transformation")
    modelcbto.setParam('OutputFlag', False)
    DictVar = {crit: modelcbto.addVar(vtype=GRB.INTEGER, lb=1) for crit in range(len(W_array))}

    for k in range(1, len(L) // 2 + 2):
        modelcbto.addConstr(
            quicksum([DictVar[crit] for crit in L[k]]) >= quicksum([DictVar[crit] for crit in L[k - 1]]) + 1)

    modelcbto.setObjective(quicksum(DictVar.values()), GRB.MINIMIZE)
    modelcbto.optimize()
    if modelcbto.status == GRB.OPTIMAL:
        # DictVal = {crit: DictVar[crit].x for crit in DictVar}
        W_array_transformed = np.zeros(len(W_array))
        for crit in DictVar:
            W_array_transformed[crit] = DictVar[crit].x
        return W_array_transformed
    else:
        raise Exception("impossible de tranformer en CBTO")


if __name__ == "__main__":

    id_save_instances = 1
    id_iteration = 0
    Inst_List = [(7, 12, 999)]      # À MODIFIER : Chaque triplet est composé de la valeur m (nbre de critères), de la
                                     # taille n de l'ensemble A d'alternatives et de la valeur d'une graine aléatoire
                                     # Attention, pour obtenir les mêmes résultats que ceux du manuscrit, conserver n = 12 et
                                     # la graine à 999.
    k_best_value = 3                 # à ne modifier que si on souhaite réaliser d'autres tests (autres que ceux du manuscrit)

    for elm in Inst_List:
        D = 0
        T = 0

        id_save_instances = 1
        id_iteration = 0
        m, xp_size_A, xp_seed = elm
        assert xp_size_A == 12
        DATA_NB_DECISIONS = [0] * 13
        DATA_NB_CHOIX_XPLAINED = [0] * 13
        DATA_NB_ELIMINATIONS_XPLAINED = [0] * 13
        np.random.seed(xp_seed)
        # for xp_size_AR in range(ceil(m / 2), floor(3 * m / 2) + 1):
        nb_conclusions_robustes = 0
        print(datetime.datetime.now())
        while nb_conclusions_robustes < 1000:

            for xp_size_AR in range(m-2, m+3):
                # print("|AR|", "=", xp_size_AR)
                omega_score, setA_as_binary_np_arrays, setAR_as_binary_np_arrays = generate_instance(m,
                                                                                                     size_A=xp_size_A,
                                                                                                     size_AR=xp_size_AR)

                # modifier omega score : remplacer avec un CBTO si m = 6
                if m == 6:
                    id_fichier = np.random.randint(1, len(os.listdir(directory)) + 1)
                    with open(directory + f'/model{id_fichier}.csv') as utilityFile:
                        reader = csv.DictReader(utilityFile)
                        w_list = list()
                        for row in reader:
                            for criterion in reader.fieldnames:
                                w_list.append(int(row[criterion]))  # rester sur KR-CBTO6

                    np.random.shuffle(w_list)

                    omega_score = np.array(w_list)
                    setA_as_binary_np_arrays = list(
                        sorted(setA_as_binary_np_arrays, key=lambda arr: sum(arr * omega_score), reverse=True))
                    setAR_as_binary_np_arrays = list(
                        sorted(setAR_as_binary_np_arrays, key=lambda arr: sum(arr * omega_score), reverse=True))
                    print("deh")
                else:
                    omega_score = transform_omega_function_into_cbto(omega_score)
                    # print(omega_score)

                PI_cofactors = [setAR_as_binary_np_arrays[i] - setAR_as_binary_np_arrays[i + 1] for i in
                                range(0, len(setAR_as_binary_np_arrays) - 1)]
                PI_statements = [(setAR_as_binary_np_arrays[i], setAR_as_binary_np_arrays[i + 1]) for i in
                                 range(0, len(setAR_as_binary_np_arrays) - 1)]

                The_statements_deduits_dict = {(k1, k2): (setA_as_binary_np_arrays[k1], setA_as_binary_np_arrays[k2])
                                               for k1 in range(len(setA_as_binary_np_arrays)) for k2 in
                                               range(k1 + 1, len(setA_as_binary_np_arrays)) if
                                               np_computation(
                                                   setA_as_binary_np_arrays[k1] - setA_as_binary_np_arrays[k2],
                                                   PI_cofactors)}

                The_explanations_dict = compute_robust_explanations(The_statements_deduits_dict, PI_statements,
                                                                    PI_cofactors)

                Schemes_used = set()
                nb_decision_choix = 0
                nb_choix_explained = 0
                LesChoix = list()
                LesChoixExpliques = list()
                LesNonChoix = list()
                LesNonChoixExpliques = list()
                # Decisions Choix
                for i_best in range(k_best_value):
                    k_best_alternative_array = setA_as_binary_np_arrays[i_best]
                    nb_necessary_dominated = 0
                    # decided = True
                    for i_worst in range(i_best + 1, len(setA_as_binary_np_arrays)):
                        k_worst_alternative_array = setA_as_binary_np_arrays[i_worst]
                        if np_computation(k_best_alternative_array - k_worst_alternative_array, PI_cofactors):
                            nb_necessary_dominated += 1
                    if nb_necessary_dominated >= xp_size_A - k_best_value:
                        nb_decision_choix += 1
                        LesChoix.append(i_best)
                        nb_necessary_dominated_explained = 0
                        TmpSchemes = set()
                        for i_worst in range(i_best + 1, len(setA_as_binary_np_arrays)):
                            if (i_best, i_worst) in The_explanations_dict:
                                nb_necessary_dominated_explained += 1
                                TmpSchemes.add(The_explanations_dict[(i_best, i_worst)][0])
                        if nb_necessary_dominated_explained >= xp_size_A - k_best_value:
                            LesChoixExpliques.append(i_best)
                            nb_choix_explained += 1
                            Schemes_used = Schemes_used | TmpSchemes

                nb_decision_elimination = 0
                nb_elimination_explained = 0
                # Decisions non choix
                for i_worst in range(k_best_value, len(setA_as_binary_np_arrays)):
                    k_worst_alternative_array = setA_as_binary_np_arrays[i_worst]
                    nb_necessary_dominating = 0
                    for i_best in range(i_worst):
                        k_best_alternative_array = setA_as_binary_np_arrays[i_best]
                        # for k_best_alternative_array in setA_as_binary_np_arrays[:k_best_value]:
                        if np_computation(k_best_alternative_array - k_worst_alternative_array, PI_cofactors):
                            nb_necessary_dominating += 1
                    if nb_necessary_dominating >= k_best_value:
                        nb_decision_elimination += 1
                        LesNonChoix.append(i_worst)
                        nb_necessary_dominating_explained = 0
                        TmpSchemes = set()
                        for i_best in range(i_worst):
                            if (i_best, i_worst) in The_explanations_dict:
                                nb_necessary_dominating_explained += 1
                                TmpSchemes.add(The_explanations_dict[(i_best, i_worst)][0])
                        if nb_necessary_dominating_explained >= k_best_value:
                            nb_elimination_explained += 1
                            LesNonChoixExpliques.append(i_worst)
                            Schemes_used = Schemes_used | TmpSchemes

                nb_conclusions_robustes += nb_decision_choix + nb_decision_elimination

                DATA_NB_DECISIONS[nb_decision_choix + nb_decision_elimination] += 1
                DATA_NB_CHOIX_XPLAINED[nb_decision_choix + nb_decision_elimination] += nb_choix_explained
                DATA_NB_ELIMINATIONS_XPLAINED[nb_decision_choix + nb_decision_elimination] += nb_elimination_explained
                print(id_save_instances, "Choix (explained)", nb_decision_choix, "(", nb_choix_explained, ")",
                      "Elimination (explained)", nb_decision_elimination, "(", nb_elimination_explained, ")",
                      "Total (explained)",
                      nb_decision_choix + nb_decision_elimination, "(", nb_choix_explained + nb_elimination_explained,
                      ")", Schemes_used)
                save_instance_pfia(
                    XP_DIR + f'm{m}/',
                    f'instance_{id_save_instances}_Niv0_choix{len(LesChoixExpliques)}_sur{len(LesChoix)}_nonchoix{len(LesNonChoixExpliques)}_sur{len(LesNonChoix)}_totaldecisions{len(LesChoixExpliques) + len(LesNonChoixExpliques)}_sur{len(LesChoix) + len(LesNonChoix)}_AR{xp_size_AR}.chapter3',
                    omega_score,
                    setA_as_binary_np_arrays, setAR_as_binary_np_arrays,
                    (The_statements_deduits_dict, LesChoix, LesChoixExpliques, LesNonChoix, LesNonChoixExpliques,
                     The_explanations_dict))
                id_save_instances += 1
                if nb_conclusions_robustes >= 1000:
                    break

            print("=====> ", sum(DATA_NB_CHOIX_XPLAINED) + sum(DATA_NB_ELIMINATIONS_XPLAINED), "explained sur", nb_conclusions_robustes, "deduites", datetime.datetime.now())
        print("Nombre d'instances en fonction du nombre de décisions\n", DATA_NB_DECISIONS)
        DENOMINATEUR = [j * DATA_NB_DECISIONS[j] for j in range(len(DATA_NB_DECISIONS))]
        ARRAY_CHOIX_ELIMINATIONS_XPLAINED = np.array(DATA_NB_CHOIX_XPLAINED) + np.array(DATA_NB_ELIMINATIONS_XPLAINED)
        PROPORTIONS = [round(100. * ARRAY_CHOIX_ELIMINATIONS_XPLAINED[j] / DENOMINATEUR[j], 1) for j in
                       range(len(DENOMINATEUR)) if DENOMINATEUR[j] != 0]
        print("Denominateur\n", DENOMINATEUR)
        print("XP\n", ARRAY_CHOIX_ELIMINATIONS_XPLAINED)
        print("Proportions\n", PROPORTIONS)
        print("=====> ", sum(DATA_NB_CHOIX_XPLAINED) + sum(DATA_NB_ELIMINATIONS_XPLAINED), "explained sur", nb_conclusions_robustes, "deduites", datetime.datetime.now())
