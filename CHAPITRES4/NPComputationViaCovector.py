from gurobipy import *

import numpy as np
from InstanceGenerator import generate_instance
from Utils import predicat_in_for_arrays
from Tools import CONSTRAINTSFEASIBILITYTOL, INTEGERFEASIBILITYTOL, EPSILON

import itertools as iter


def np_computation(cofactor_to_check, PI_cofactors_List):
    """cofactor_to_check : cofacteur (de type array) de la paire (x, y) dont on veut verifier l'appartenance à la N_PI
       PI_cofactors_List : Liste des cofacteurs de la PI (rangement de l'ensemble AR)"""
    model = Model("Necessary Preference Computation via Cofactors")
    model.setParam('OutputFlag', False)
    model.Params.IntFeasTol = INTEGERFEASIBILITYTOL
    model.Params.FeasibilityTol = CONSTRAINTSFEASIBILITYTOL

    # VarM : Une variable par critère pour l'utilisation de la dominance
    nb_criteria_considered_alias_m = len(cofactor_to_check)
    VarM = np.array([model.addVar(vtype=GRB.CONTINUOUS, lb=0, name="M_{}".format(i)) for i in
                     range(nb_criteria_considered_alias_m)])

    # VarL : une variable par element de la PI
    VarL = [model.addVar(vtype=GRB.CONTINUOUS, lb=0, name="L_pi_{}".format(k)) for k in range(len(PI_cofactors_List))]

    # VarLDup : instance de chaque variable dans VarL dupliquée autant de fois qu'il y a de critères afin d'effectuer la multiplication avec le cofacteur
    VarLDup = [np.array([vl] * nb_criteria_considered_alias_m) for vl in VarL]
    VarLDupMultEachPI_cofactor = [VarLDup[k] * PI_cofactors_List[k] for k in range(len(PI_cofactors_List))]

    RightMember_1 = [quicksum([elm[i] for elm in VarLDupMultEachPI_cofactor]) for i in
                     range(nb_criteria_considered_alias_m)]
    RightMember_2 = VarM

    model.update()
    LeftMember = cofactor_to_check
    for i in range(nb_criteria_considered_alias_m):
        model.addConstr(RightMember_1[i] + RightMember_2[i] == LeftMember[i], name=f'Constraint_{i}')
        # model.addConstr(LeftMember[i] == RightMember_1[i] + RightMember_2[i], name=f'Constraint_{i}') # génère une error : gurobi ne supporte pas que la constante soit mise à gauche

    model.update()
    model.optimize()

    return model.status == GRB.OPTIMAL


def is_2Order_necessary_swap_explainable(cofactor_to_check, PI_cofactors_List):
    nb_criteria_considered_alias_m = len(cofactor_to_check)
    SigmaP = [i for i in range(nb_criteria_considered_alias_m) if cofactor_to_check[i] == 1]
    SigmaD = [i for i in range(nb_criteria_considered_alias_m) if cofactor_to_check[i] == -1]

    edgeCoeffDict = dict()
    concernedBooleanVarDict = {i: list() for i in set(SigmaD) | set(SigmaP)}
    booleanVarDict = dict()

    for i in SigmaP:
        for j in SigmaD:
            fictious_cofactor = np.zeros(nb_criteria_considered_alias_m)
            fictious_cofactor[i] = 1
            fictious_cofactor[j] = -1
            if np_computation(fictious_cofactor, PI_cofactors_List):
                edgeCoeffDict[(i, j)] = 1
            else:
                edgeCoeffDict[(i, j)] = 0

    matching_gurobi_model = Model("Explain Matching Model")
    matching_gurobi_model.setParam('OutputFlag', False)
    matching_gurobi_model.Params.IntFeasTol = INTEGERFEASIBILITYTOL
    for i, j in edgeCoeffDict:
        var_ij = matching_gurobi_model.addVar(vtype=GRB.BINARY, name="B_{}_{}".format(i, j))
        concernedBooleanVarDict[i].append(var_ij)
        concernedBooleanVarDict[j].append(var_ij)
        booleanVarDict[(i, j)] = var_ij

    matching_gurobi_model.update()
    for i, VarList in concernedBooleanVarDict.items():
        matching_gurobi_model.addConstr(quicksum(VarList) <= 1)

    matching_gurobi_model.setObjective(
        quicksum([edgeCoeffDict[(i, j)] * booleanVarDict[(i, j)] for i, j in edgeCoeffDict]), GRB.MAXIMIZE)
    matching_gurobi_model.update()

    matching_gurobi_model.optimize()
    explainable = round(matching_gurobi_model.objVal) == len(SigmaD)
    return explainable


def is_2Order_mixed_swap_explainable(cofactor_to_check, PI_cofactors_List, Reco_cofactors_List):
    nb_criteria_considered_alias_m = len(cofactor_to_check)
    model = Model("Mixed Explanations DA2PL20")
    model.setParam('OutputFlag', False)
    model.Params.FeasibilityTol = CONSTRAINTSFEASIBILITYTOL
    model.Params.IntFeasTol = INTEGERFEASIBILITYTOL
    Var_omega = [model.addVar(vtype=GRB.CONTINUOUS, lb=0, name="w_{}".format(i)) for i in
                 range(nb_criteria_considered_alias_m)]

    pros = {i for i in range(nb_criteria_considered_alias_m) if cofactor_to_check[i] == 1}
    cons = {i for i in range(nb_criteria_considered_alias_m) if cofactor_to_check[i] == -1}

    Var_Bij = {(i, j): model.addVar(vtype=GRB.BINARY, name=f'b_{i}_{j}') for i in pros for j in cons}

    PI_and_Reco_cofactors_List = PI_cofactors_List + Reco_cofactors_List
    # Constraint (6) in DA2PL20
    for elmt in PI_and_Reco_cofactors_List:
        model.addConstr(quicksum(np.array(Var_omega) * elmt) >= 0)
        # model.addConstr(quicksum(np.array(Var_omega) * elmt) >= CONSTRAINTSFEASIBILITYTOL)

    # Constraint (7) in DA2PL20
    model.addConstr(quicksum(Var_omega) == 1.)

    # Constraint (9) in DA2PL20
    for j in cons:
        model.addConstr(quicksum([Var_Bij[(i_, j)] for i_ in pros]) == 1)

    # Constraint (10) in DA2PL20
    for i in pros:
        model.addConstr(quicksum([Var_Bij[(i, j_)] for j_ in cons]) <= 1)

    # Constraint (11) in DA2PL20
    for (i, j), var_bij in Var_Bij.items():
        model.addConstr(Var_omega[i] - Var_omega[j] >= var_bij - 1)
        # model.addConstr(Var_omega[i] - Var_omega[j] >= (1 + CONSTRAINTSFEASIBILITYTOL) * var_bij - 1)

    objFunction = LinExpr()
    # Liste_swap_necessaire = list()
    for i in pros:
        for j in cons:
            fictious_cofactor = np.zeros(nb_criteria_considered_alias_m)
            fictious_cofactor[i] = 1
            fictious_cofactor[j] = -1
            if np_computation(fictious_cofactor, PI_cofactors_List):
                objFunction += Var_Bij[(i, j)]
                # Liste_swap_necessaire.append((i, j))

    model.update()
    model.setObjective(objFunction, GRB.MAXIMIZE)
    model.optimize()

    if model.status == GRB.OPTIMAL:
        # construire cofactor possible
        all_swaps_cofactor_Liste = list()
        for (i, j), var_bij in Var_Bij.items():
            # if (i, j) not in Liste_swap_necessaire and round(var_bij.x) == 1:
            if round(var_bij.x) == 1:
                fictious_cofactor = np.zeros(nb_criteria_considered_alias_m)
                fictious_cofactor[i] = 1
                fictious_cofactor[j] = -1
                all_swaps_cofactor_Liste.append(fictious_cofactor)
        # ICI, IL FAUT QUE lA LISTE SOIT NON VIDE; SI VIDE PEUT-ÊTRE RETOURNER FALSE
        if len(all_swaps_cofactor_Liste) == 0:
            return False, None, list()
            # sys.stdout.write(str([Var_omega[i].x for i in range(nb_criteria_considered_alias_m)]))
            # sys.stdout.write(str({(i, j): var_bij.x for (i, j), var_bij in Var_Bij.items()}))
            # sys.stdout.write(f'pros {str(pros)}, cons  {str(cons)}')
            # sys.stdout.write(f'Necessaire, {Liste_swap_necessaire}')
        return True, len(cons) - round(model.objVal), all_swaps_cofactor_Liste
    return False, None, list()


def is_2Order_mixed_swap_explainable_HEURISTICS(cofactor_to_check, PI_cofactors_List, Reco_cofactors_List):
    nb_criteria_considered_alias_m = len(cofactor_to_check)
    pros = {i for i in range(nb_criteria_considered_alias_m) if cofactor_to_check[i] == 1}
    cons = {i for i in range(nb_criteria_considered_alias_m) if cofactor_to_check[i] == -1}
    All_pairs_pros_cons = [(i, j) for i in pros for j in cons]
    Dico_combi_cofactors_list = dict()
    for combi in iter.combinations(All_pairs_pros_cons, len(cons)):
        if len(set([elmt[1] for elmt in combi])) == len(cons) and len(set([elmt[0] for elmt in combi])) == len(
                [elmt[0] for elmt in combi]):
            corresponding_cofactors = list()
            for (i, j) in combi:
                fictious_cofactor = np.zeros(nb_criteria_considered_alias_m)
                fictious_cofactor[i] = 1
                fictious_cofactor[j] = -1
                corresponding_cofactors.append(fictious_cofactor)
            Dico_combi_cofactors_list[combi] = (
                proportion_possible(corresponding_cofactors, PI_cofactors_List + Reco_cofactors_List,
                                    nb_tirs_in_PI_RECO=100), corresponding_cofactors)

    better_combi_key = max(Dico_combi_cofactors_list.keys(), key=lambda item: item[0][0])
    better_combi = Dico_combi_cofactors_list[better_combi_key]
    # print(better_combi[0])
    if better_combi[0][0] == 0:
        print("wooooooooooooooooo")
        return is_2Order_mixed_swap_explainable(cofactor_to_check, PI_cofactors_List, Reco_cofactors_List)
    else:
        nb_non_necessaires = 0
        for cofac in better_combi[1]:
            if not np_computation(cofac, PI_cofactors_List):
                nb_non_necessaires += 1
        return True, nb_non_necessaires, better_combi[1]


if __name__ == "__main__":
    nb_iterations = 1000
    cpt_deduits = 0
    cpt_dans_PI = 0
    cpt_deduits_non_PI_eligibles = 0
    cpt_deduits_non_PI_eligibles_not_trivial = 0
    cpt_deduits_non_PI_eligibles_explicables_necessary = 0
    cpt_gain_mixed = 0
    cpt_gain_mixed_exactement_1_possible = 1
    cpt_gain_mixed_exactement_2_possible = 2
    Liste_value_p = list()
    while nb_iterations > 0:
        nb_iterations -= 1
        omega_score, setA_as_binary_np_arrays, setAR_as_binary_np_arrays = generate_instance(15, size_A=10, size_AR=10)
        reference_cofactors = [setAR_as_binary_np_arrays[i] - setAR_as_binary_np_arrays[i + 1] for i in
                               range(0, len(setAR_as_binary_np_arrays) - 1)]
        reco_cofactors = [setA_as_binary_np_arrays[0] - setA_as_binary_np_arrays[k] for k in
                          range(1, len(setA_as_binary_np_arrays))]
        for challenger in setA_as_binary_np_arrays[1:]:  # CHOICE
            inst_cofactor_to_check = setA_as_binary_np_arrays[0] - challenger
            # for k in range(len(setA_as_binary_np_arrays) - 1):                   # RANKING
            #     inst_cofactor_to_check = setA_as_binary_np_arrays[k] - setA_as_binary_np_arrays[k+1]
            if predicat_in_for_arrays(setA_as_binary_np_arrays[0],  # CHOICE
                                      setAR_as_binary_np_arrays) and predicat_in_for_arrays(challenger,
                                                                                            setAR_as_binary_np_arrays):
                # if predicat_in_for_arrays(setA_as_binary_np_arrays[k],                             # RANKING
                #                           setAR_as_binary_np_arrays) and predicat_in_for_arrays(setA_as_binary_np_arrays[k+1],
                #                                                                                 setAR_as_binary_np_arrays):
                cpt_deduits += 1
                cpt_dans_PI += 1
            elif np_computation(inst_cofactor_to_check, reference_cofactors):
                cpt_deduits += 1
                nb_pros = np.count_nonzero(inst_cofactor_to_check == 1)
                nb_cons = np.count_nonzero(inst_cofactor_to_check == -1)
                if nb_pros >= nb_cons:  # ELIGIBLE
                    cpt_deduits_non_PI_eligibles += 1
                    if nb_pros >= 2:
                        cpt_deduits_non_PI_eligibles_not_trivial += 1
                    if is_2Order_necessary_swap_explainable(inst_cofactor_to_check, reference_cofactors):
                        cpt_deduits_non_PI_eligibles_explicables_necessary += 1
                    else:
                        # res, nb_possible_swap, All_Swaps_Cofactor_List = is_2Order_mixed_swap_explainable(
                        #     inst_cofactor_to_check,
                        #     reference_cofactors, reco_cofactors)

                        res, nb_possible_swap, All_Swaps_Cofactor_List = is_2Order_mixed_swap_explainable_HEURISTICS(
                            inst_cofactor_to_check,
                            reference_cofactors, reco_cofactors)

                        # if nb_possible_swap != len(Possible_Cofactor_List):
                        #     print(nb_possible_swap, Possible_Cofactor_List, len(Possible_Cofactor_List))
                        # assert nb_possible_swap == len(Possible_Cofactor_List)         # CHECKER POURQUOI QUELQUES FOIS ÇA CRACHE

                        if res:
                            cpt_gain_mixed += 1
                            value_p, total_considered, nb_attempts = proportion_possible(All_Swaps_Cofactor_List,
                                                                                         reference_cofactors)
                            Liste_value_p.append(value_p)
                            if nb_possible_swap == 1:
                                cpt_gain_mixed_exactement_1_possible += 1
                            if nb_possible_swap == 2:
                                cpt_gain_mixed_exactement_2_possible += 1
                            if nb_possible_swap == 0:
                                print("CURIEUX")

        if (nb_iterations + 1) % 20 == 0:
            print(nb_iterations, "cpt", cpt_deduits)
    print("Compteur :", cpt_deduits, " dont ", cpt_dans_PI, " dans la PI")
    print("Déduits non PI eligibles", cpt_deduits_non_PI_eligibles, " dont ", cpt_deduits_non_PI_eligibles_not_trivial,
          "non triviaux")
    print("Déduits non PI eligibles", cpt_deduits_non_PI_eligibles, " dont ",
          cpt_deduits_non_PI_eligibles_explicables_necessary,
          "explicables necessairement",
          round(100 * cpt_deduits_non_PI_eligibles_explicables_necessary / cpt_deduits_non_PI_eligibles, 2), "%")
    print("Gain possible avec 1 swap possible", "valcpt", cpt_gain_mixed_exactement_1_possible, "pource",
          round(100 * cpt_gain_mixed_exactement_1_possible / cpt_deduits_non_PI_eligibles, 2), "%")
    print("Gain possible avec 2 swap possible", "valcpt", cpt_gain_mixed_exactement_2_possible, "pource",
          round(100 * cpt_gain_mixed_exactement_2_possible / cpt_deduits_non_PI_eligibles, 2), "%")

    print("Gain possible general", "valcpt", cpt_gain_mixed, "pource",
          round(100 * cpt_gain_mixed / cpt_deduits_non_PI_eligibles, 2), "%")

    print("Hypervolume", "Min", min(Liste_value_p), "max", max(Liste_value_p), "moy", np.mean(Liste_value_p))
