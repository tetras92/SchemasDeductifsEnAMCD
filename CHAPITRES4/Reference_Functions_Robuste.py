import numpy as np
from gurobipy import *
from Tools import INTEGERFEASIBILITYTOL
from NPComputationViaCovector import np_computation
import itertools as itert
import datetime

PFIA_XP_DIR = f'/home/manuel239/PycharmProjects/MOMA/CORE/IJAR/XP/PFIA2023/XP1Best/'


def doPa_robust(covector_to_check, step=0, BDB=None):
    if BDB is None:
        BDB = dict()
    if tuple(covector_to_check) in BDB['[DoPa]']:
        return BDB['[DoPa]'][tuple(covector_to_check)]

    prefix = "" if step == 0 else "->"
    if all(covector_to_check >= 0) and not (all(covector_to_check == 0)):
        BDB['[DoPa]'][tuple(covector_to_check)] = (True, (prefix + "[DoPa]", None, None, None))
        return True, (prefix + "[DoPa]", None, None, None)

    BDB['[DoPa]'][tuple(covector_to_check)] = False, (None, None, None, None)
    return False, (None, None, None, None)


def cePa_robust(to_explain_statement, PI_statements_List, PI_cofactors_List, step=0, BDB=None):
    if BDB is None:
        BDB = dict()

    prefix = "" if step == 0 else "->"
    alt_x, alt_y = to_explain_statement[0], to_explain_statement[1]
    if (tuple(alt_x), tuple(alt_y)) in BDB['[CePa]']:
        return BDB['[CePa]'][(tuple(alt_x), tuple(alt_y))]

    cofactor_to_check = alt_x - alt_y
    for alt_u, alt_v in PI_statements_List:
        if all(alt_u - alt_x == 0) and all(alt_v - alt_y == 0):
            BDB['[CePa]'][(tuple(alt_x), tuple(alt_y))] = False, (None, None, None, None)
            return False, (None, None, None, None)  # CePa exclut l'identite : 13/07/2023

    for k in range(len(PI_cofactors_List)):
        cof_pi = PI_cofactors_List[k]
        if all(cofactor_to_check - cof_pi == 0):
            details = (prefix + "[CePa]", ("".join(
                [chr(ord('a') + i) for i in range(len(PI_statements_List[k][0])) if PI_statements_List[k][0][i] == 1]),
                                              "".join(
                                                  [chr(ord('a') + i) for i in range(len(PI_statements_List[k][1])) if
                                                   PI_statements_List[k][1][i] == 1])), None, None)
            BDB['[CePa]'][(tuple(alt_x), tuple(alt_y))] = True, details
            return True, details

    BDB['[CePa]'][(tuple(alt_x), tuple(alt_y))] = False, (None, None, None, None)
    return False, (None, None, None, None)


def tr_robust(to_explain_statement, PI_statements_List, step=0, BDB=None):
    if BDB is None:
        BDB = dict()

    prefix = "" if step == 0 else "->"
    alt_x, alt_y = to_explain_statement[0], to_explain_statement[1]
    if (tuple(alt_x), tuple(alt_y)) in BDB['[Tr]']:
        return BDB['[Tr]'][(tuple(alt_x), tuple(alt_y))]

    Successeur_Dict = {st_k: [st_k_prime for st_k_prime in range(len(PI_statements_List)) if
                              np.all(PI_statements_List[st_k][1] == PI_statements_List[st_k_prime][0])] + (
                                 [-2] if np.all(PI_statements_List[st_k][1] == alt_y) else list()) for st_k in
                       range(len(PI_statements_List))}
    Successeur_Dict[-1] = [st_k for st_k in range(len(PI_statements_List)) if
                           np.all(PI_statements_List[st_k][0] == alt_x)]
    Successeur_Dict[-2] = list()

    model_gurobi = Model("Calcul PCC for transitive explanations")
    model_gurobi.setParam('OutputFlag', False)

    Arcs_Liaisons_Var_Dict = {(elmt, succ_elmt): model_gurobi.addVar(vtype=GRB.BINARY, name=f'b_{elmt}_{succ_elmt}') for
                              elmt, LSE in Successeur_Dict.items() for succ_elmt in LSE}

    # contrainte : exactement un arc sortant du sommet -1
    model_gurobi.addConstr(quicksum(
        [var_couple_pred_succ for couple_pred_succ, var_couple_pred_succ in Arcs_Liaisons_Var_Dict.items() if
         couple_pred_succ[0] == -1]) == 1)

    # contraintes : continuité du chemin
    for st_k in range(len(PI_statements_List)):
        model_gurobi.addConstr(
            quicksum(
                [var_couple_pred_succ for couple_pred_succ, var_couple_pred_succ in Arcs_Liaisons_Var_Dict.items() if
                 couple_pred_succ[0] == st_k]) == quicksum(
                [var_couple_pred_succ for couple_pred_succ, var_couple_pred_succ in Arcs_Liaisons_Var_Dict.items() if
                 couple_pred_succ[1] == st_k]))

    model_gurobi.update()

    model_gurobi.setObjective(quicksum(Arcs_Liaisons_Var_Dict.values()), GRB.MINIMIZE)
    model_gurobi.optimize()
    if model_gurobi.status == GRB.OPTIMAL:
        Chaine_des_st_k = [-2]
        elmt_courant = Chaine_des_st_k[-1]
        while elmt_courant != -1:
            for couple_pred_succ, var_couple_pred_succ in Arcs_Liaisons_Var_Dict.items():
                if couple_pred_succ[1] == elmt_courant and round(var_couple_pred_succ.x) == 1:
                    Chaine_des_st_k.append(couple_pred_succ[0])
                    break
            elmt_courant = Chaine_des_st_k[-1]
        Chaine_des_st_k.reverse()
        details = (prefix + "[Tr]",
                   Chaine_des_st_k[1:-1],
                   [("".join(
                       [chr(ord('a') + i) for i in range(len(PI_statements_List[a][0])) if
                        PI_statements_List[a][0][i] == 1]),
                     "".join(
                         [chr(ord('a') + i) for i in range(len(PI_statements_List[a][1])) if
                          PI_statements_List[a][1][i] == 1]))
                       for a
                       in Chaine_des_st_k[1:-1]],
                   len(Chaine_des_st_k[1:-1]))
        BDB['[Tr]'][(tuple(alt_x), tuple(alt_y))] = (True, details)
        return True, details
    BDB['[Tr]'][(tuple(alt_x), tuple(alt_y))] = False, (None, None, None, None)
    return False, (None, None, None, None)


def trE_robust(to_explain_statement, PI_statements_List, step=0, BDB=None):
    if BDB is None:
        BDB = dict()

    prefix = "" if step == 0 else "->"
    alt_x, alt_y = to_explain_statement[0], to_explain_statement[1]
    if (tuple(alt_x), tuple(alt_y)) in BDB['[TrE]']:
        return BDB['[TrE]'][(tuple(alt_x), tuple(alt_y))]

    def pareto_domine(dominant_array, dominated_array):
        return all(dominant_array - dominated_array >= 0) and not (
            all(dominant_array - dominated_array == 0))  # maj ce 14/06/2023

    Successeur_Dict = {st_k: [st_k_prime for st_k_prime in range(len(PI_statements_List)) if
                              np.all(PI_statements_List[st_k][1] == PI_statements_List[st_k_prime][0])
                              or pareto_domine(PI_statements_List[st_k][1], PI_statements_List[st_k_prime][0])] + (
                                 [-2] if pareto_domine(PI_statements_List[st_k][1], alt_y) or np.all(
                                     PI_statements_List[st_k][1] == alt_y) else list()) for st_k in
                       range(len(PI_statements_List))}
    Successeur_Dict[-1] = [st_k for st_k in range(len(PI_statements_List)) if
                           pareto_domine(alt_x, PI_statements_List[st_k][0]) or np.all(
                               alt_x == PI_statements_List[st_k][0])]
    # + (
    #     [-2] if pareto_domine(alt_x, alt_y) or np.all(alt_x == alt_y) else [])  # ça c'est du dopa ou id
    Successeur_Dict[-2] = list()

    model_gurobi = Model("Calcul PCC for transitive Pareto augmented explanations")
    model_gurobi.setParam('OutputFlag', False)

    Arcs_Liaisons_Var_Dict = {(elmt, succ_elmt): model_gurobi.addVar(vtype=GRB.BINARY, name=f'b_{elmt}_{succ_elmt}') for
                              elmt, LSE in Successeur_Dict.items() for succ_elmt in LSE}

    # contrainte : exactement un arc sortant du sommet -1
    model_gurobi.addConstr(quicksum(
        [var_couple_pred_succ for couple_pred_succ, var_couple_pred_succ in Arcs_Liaisons_Var_Dict.items() if
         couple_pred_succ[0] == -1]) == 1)

    # contraintes : continuité du chemin
    for st_k in range(len(PI_statements_List)):
        model_gurobi.addConstr(
            quicksum(
                [var_couple_pred_succ for couple_pred_succ, var_couple_pred_succ in Arcs_Liaisons_Var_Dict.items() if
                 couple_pred_succ[0] == st_k]) == quicksum(
                [var_couple_pred_succ for couple_pred_succ, var_couple_pred_succ in Arcs_Liaisons_Var_Dict.items() if
                 couple_pred_succ[1] == st_k]))

    model_gurobi.update()

    model_gurobi.setObjective(quicksum(Arcs_Liaisons_Var_Dict.values()), GRB.MINIMIZE)
    model_gurobi.optimize()
    if model_gurobi.status == GRB.OPTIMAL:
        Chaine_des_st_k = [-2]
        elmt_courant = Chaine_des_st_k[-1]
        while elmt_courant != -1:
            for couple_pred_succ, var_couple_pred_succ in Arcs_Liaisons_Var_Dict.items():
                if couple_pred_succ[1] == elmt_courant and round(var_couple_pred_succ.x) == 1:
                    Chaine_des_st_k.append(couple_pred_succ[0])
                    break
            elmt_courant = Chaine_des_st_k[-1]
        # print("iiiiii")
        Chaine_des_st_k.reverse()
        PresenceDominance = [None]
        for i_pos in range(1, len(Chaine_des_st_k)):
            if Chaine_des_st_k[i_pos] != -2:
                statement_pos = PI_statements_List[Chaine_des_st_k[i_pos]]
                dominant, dominated = statement_pos
                if Chaine_des_st_k[i_pos - 1] == -1:
                    to_consider = alt_x
                else:
                    # to_consider = PI_statements_List[i_pos - 1][1] BUG CORRIGE CE 16 AOUT
                    to_consider = PI_statements_List[Chaine_des_st_k[i_pos - 1]][1]
                if pareto_domine(to_consider, dominant):
                    PresenceDominance.append("D")
                else:
                    PresenceDominance.append(None)
            else:
                dominated = alt_y
                if Chaine_des_st_k[i_pos - 1] == -1:
                    to_consider = alt_x
                else:
                    to_consider = PI_statements_List[Chaine_des_st_k[i_pos - 1]][1]
                if pareto_domine(to_consider, dominated):
                    PresenceDominance.append("D")
                else:
                    PresenceDominance.append(None)
        Chaine_a_retourner = list()
        assert len(Chaine_des_st_k) == len(PresenceDominance)
        for i_pos in range(1, len(Chaine_des_st_k)):
            if i_pos != len(Chaine_des_st_k) - 1:
                if PresenceDominance[i_pos] == "D":
                    Chaine_a_retourner.append("D")
                Chaine_a_retourner.append(("".join(
                    [chr(ord('a') + i) for i in range(len(PI_statements_List[Chaine_des_st_k[i_pos]][0])) if
                     PI_statements_List[Chaine_des_st_k[i_pos]][0][i] == 1]),
                                           "".join([chr(ord('a') + i) for i in
                                                    range(len(PI_statements_List[Chaine_des_st_k[i_pos]][1])) if
                                                    PI_statements_List[Chaine_des_st_k[i_pos]][1][i] == 1])))
            elif PresenceDominance[i_pos] == "D":
                Chaine_a_retourner.append("D")
        details = (prefix + "[TrE]",
                   Chaine_des_st_k[1:-1],
                   Chaine_a_retourner,
                   len(Chaine_des_st_k[1:-1]))
        BDB['[TrE]'][(tuple(alt_x), tuple(alt_y))] = True, details
        return True, details

    BDB['[TrE]'][(tuple(alt_x), tuple(alt_y))] = (False, (None, None, None, None))
    return False, (None, None, None, None)


def cov0_robust(covector_to_check, PI_cofactors_List, step=0, BDB=None):
    if BDB is None:
        BDB = dict()
    # print(covector_to_check)
    # if tuple(covector_to_check) in BDB['[Cov]']:
    #     return BDB['[Cov]'][tuple(covector_to_check)]

    prefix = "" if step == 0 else "->"
    nb_criteria_considered_alias_m = len(covector_to_check)
    SigmaP = [i_ for i_ in range(nb_criteria_considered_alias_m) if covector_to_check[i_] == 1]
    SigmaD = [i_ for i_ in range(nb_criteria_considered_alias_m) if covector_to_check[i_] == -1]

    edgeCoeffDict = dict()
    concernedBooleanVarDict = {i: list() for i in set(SigmaD) | set(SigmaP)}
    booleanVarDict = dict()

    for i_ in SigmaP:
        for j in SigmaD:
            fictious_cofactor = np.zeros(nb_criteria_considered_alias_m)
            fictious_cofactor[i_] = 1
            fictious_cofactor[j] = -1
            if np_computation(fictious_cofactor, PI_cofactors_List):
                edgeCoeffDict[(i_, j)] = 1
            else:
                edgeCoeffDict[(i_, j)] = 0

    matching_gurobi_model = Model("Explain Matching Model")
    matching_gurobi_model.setParam('OutputFlag', False)
    matching_gurobi_model.Params.IntFeasTol = INTEGERFEASIBILITYTOL
    for i_, j in edgeCoeffDict:
        var_ij = matching_gurobi_model.addVar(vtype=GRB.BINARY, name="B_{}_{}".format(i_, j))
        concernedBooleanVarDict[i_].append(var_ij)
        concernedBooleanVarDict[j].append(var_ij)
        booleanVarDict[(i_, j)] = var_ij

    matching_gurobi_model.update()
    for i_, VarList in concernedBooleanVarDict.items():
        matching_gurobi_model.addConstr(quicksum(VarList) <= 1)

    matching_gurobi_model.setObjective(
        quicksum([edgeCoeffDict[(i, j)] * booleanVarDict[(i, j)] for i, j in edgeCoeffDict]), GRB.MAXIMIZE)
    matching_gurobi_model.update()

    matching_gurobi_model.optimize()
    explainable = round(matching_gurobi_model.objVal) == len(SigmaD)
    if not explainable:
        # BDB['[Cov]'][tuple(covector_to_check)] = False, (None, None, None, None)
        return False, (None, None, None, None)

    details = (
        prefix + "[Cov0]", [({chr(ord('a') + i_)}, {chr(ord('a') + j)}) for (i_, j), var_ij in booleanVarDict.items() if
                            round(var_ij.x) == 1] + (
            # [] if len(SigmaP) == len(SigmaD) else [(set(SigmaP) - set(SigmaD),)]), None,                                       # bug recuperation resultat
            [] if len(SigmaP) == len(SigmaD) else [
                (set([chr(ord('a') + pr) for pr in SigmaP]) - {chr(ord('a') + i_) for (i_, j), var_ij in
                                                               booleanVarDict.items() if
                                                               round(var_ij.x) == 1},)]), None,
        # bug recuperation resultat
        len(SigmaD) + (1 if len(SigmaP) > len(SigmaD) else 0))
    # BDB['[Cov]'][tuple(covector_to_check)] = True, details

    return True, details


def covD0_robust(to_explain_statement, PI_statements_List,
                 PI_cofactors_List, step=0, BDB=None):
    if BDB is None:
        BDB = dict()

    if (tuple(to_explain_statement[0]), tuple(to_explain_statement[1])) in BDB['[CovD]']:
        return BDB['[CovD]'][(tuple(to_explain_statement[0]), tuple(to_explain_statement[1]))]

    prefix = "" if step == 0 else "->"
    # identifier les statements du milieu
    identifiant_statement_milieu_List = list()
    for st_selected_ in range(len(PI_statements_List)):
        cofactor_x_z1 = to_explain_statement[0] - PI_statements_List[st_selected_][0]
        cofactor_z2_y = PI_statements_List[st_selected_][1] - to_explain_statement[1]
        if np_computation(cofactor_x_z1, PI_cofactors_List) and np_computation(cofactor_z2_y, PI_cofactors_List):
            identifiant_statement_milieu_List.append(st_selected_)

    if len(identifiant_statement_milieu_List) == 0:
        BDB['[CovD]'][(tuple(to_explain_statement[0]), tuple(to_explain_statement[1]))] = False, (None, None, None, None)
        return False, (None, None, None, None)

    Graphe_distances = dict()  # id de x = -1 ; id de y = -2
    Explanations_dict = dict()
    for st_du_milieu in identifiant_statement_milieu_List:
        statement_courant = PI_statements_List[st_du_milieu]
        # distance x -> z1
        cofactor_x_z1_courant = to_explain_statement[0] - statement_courant[0]
        is_kb_explainable_, details_ = cov0_robust(
            cofactor_x_z1_courant, PI_cofactors_List, step=step, BDB=BDB)
        if is_kb_explainable_:
            Graphe_distances[(-1, st_du_milieu)] = details_[
                -1]
            Explanations_dict[(-1, st_du_milieu)] = details_
        # distance z2 -> y
        cofactor_z2_y_courant = statement_courant[1] - to_explain_statement[1]
        is_kb_explainable_, details_ = cov0_robust(
            cofactor_z2_y_courant, PI_cofactors_List, step=step, BDB=BDB)
        if is_kb_explainable_:
            Graphe_distances[(st_du_milieu, -2)] = details_[-1]
            Explanations_dict[(st_du_milieu, -2)] = details_

    # Completer graphe des distances avec les paires de statements du milieu
    for st_du_milieu1 in identifiant_statement_milieu_List:
        for st_du_milieu3 in identifiant_statement_milieu_List:
            if st_du_milieu1 != st_du_milieu3:
                statement1_2 = PI_statements_List[st_du_milieu1]
                statement3_4 = PI_statements_List[st_du_milieu3]
                cofactor_2_3 = statement1_2[1] - statement3_4[0]
                cofactor_4_1 = statement3_4[1] - statement1_2[0]
                is_kb_explainable_2_3, details_2_3 = cov0_robust(
                    cofactor_2_3, PI_cofactors_List, step=step, BDB=BDB)
                is_kb_explainable_4_1, details_4_1 = cov0_robust(
                    cofactor_4_1, PI_cofactors_List, step=step, BDB=BDB)
                if is_kb_explainable_2_3:
                    Graphe_distances[(st_du_milieu1, st_du_milieu3)] = details_2_3[-1]
                    Explanations_dict[(st_du_milieu1, st_du_milieu3)] = details_2_3
                if is_kb_explainable_4_1:
                    Graphe_distances[(st_du_milieu3, st_du_milieu1)] = details_4_1[-1]
                    Explanations_dict[(st_du_milieu3, st_du_milieu1)] = details_4_1
    # Calcul du PCC
    model_gurobi = Model("Calcul PCC for explanations")
    model_gurobi.setParam('OutputFlag', False)
    Arcs_Liaisons_Var_Dict = {arc: model_gurobi.addVar(vtype=GRB.BINARY, name=f'b_{arc}') for arc in Graphe_distances}

    # contrainte : exactement un arc sortant du sommet -1
    model_gurobi.addConstr(quicksum([var_b for arc, var_b in Arcs_Liaisons_Var_Dict.items() if arc[0] == -1]) == 1)

    # contraintes : continuité du chemin
    for st_du_milieu in identifiant_statement_milieu_List:
        model_gurobi.addConstr(
            quicksum([var_b for arc, var_b in Arcs_Liaisons_Var_Dict.items() if arc[0] == st_du_milieu]) == quicksum(
                [var_b for arc, var_b in Arcs_Liaisons_Var_Dict.items() if arc[1] == st_du_milieu]))

    model_gurobi.update()
    model_gurobi.setObjective(
        quicksum([(1 + Graphe_distances[arc]) * var_b for arc, var_b in Arcs_Liaisons_Var_Dict.items()]),
        GRB.MINIMIZE)
    model_gurobi.optimize()

    if model_gurobi.status == GRB.OPTIMAL:
        List_Arcs_PCC = [-1]
        courant = List_Arcs_PCC[-1]
        while courant != -2:
            suiv_courant = -2
            for st_du_milieu in identifiant_statement_milieu_List:
                if (courant, st_du_milieu) in Arcs_Liaisons_Var_Dict and round(
                        Arcs_Liaisons_Var_Dict[(courant, st_du_milieu)].x) == 1:
                    suiv_courant = st_du_milieu
                    List_Arcs_PCC.append(suiv_courant)

            courant = suiv_courant
        List_Arcs_PCC.append(-2)

        explList = list()
        for j in range(1, len(List_Arcs_PCC)):
            explList.append(Explanations_dict[(List_Arcs_PCC[j - 1], List_Arcs_PCC[j])])
            if List_Arcs_PCC[j] != -2:
                explList.append(
                    ("".join([chr(ord('a') + i) for i in range(len(PI_statements_List[List_Arcs_PCC[j]][0])) if
                              PI_statements_List[List_Arcs_PCC[j]][0][i] == 1]),
                     "".join([chr(ord('a') + i) for i in range(len(PI_statements_List[List_Arcs_PCC[j]][1])) if
                              PI_statements_List[List_Arcs_PCC[j]][1][
                                  i] == 1])))  # Mettre exactement la comparaison par paire
            else:
                explList.append(List_Arcs_PCC[j])
            # explList.append(List_Arcs_PCC[j])
        assert explList[-1] == -2

        BDB['[CovD]'][(tuple(to_explain_statement[0]), tuple(to_explain_statement[1]))] = True, (prefix + "[CovD0]", None, explList[:-1], model_gurobi.objVal - 1)
        return True, (prefix + "[CovD0]", None, explList[:-1], model_gurobi.objVal - 1)

    BDB['[CovD]'][(tuple(to_explain_statement[0]), tuple(to_explain_statement[1]))] = False, (None, None, None, None)
    return False, (None, None, None, None)


def __tr_trE_covD0_robust(to_explain_statement, PI_statements_List,
                          PI_cofactors_List, step=0, BDB=None):
    if BDB is None:
        BDB = dict()
    is_tr, details_tr = tr_robust(to_explain_statement, PI_statements_List, step=step, BDB=BDB)
    if is_tr:
        return True, details_tr
    is_tre, details_tre = trE_robust(to_explain_statement, PI_statements_List, step=step, BDB=BDB)
    if is_tre:
        return True, details_tre
    return covD0_robust(to_explain_statement, PI_statements_List,
                        PI_cofactors_List, step=step, BDB=BDB)


def un_Cg_robust(to_explain_statement, PI_statements_List, PI_cofactors_List, step=0, BDB=None):
    if BDB is None:
        BDB = dict()

    if (tuple(to_explain_statement[0]), tuple(to_explain_statement[1])) in BDB['[1-Cg]']:
        return BDB['[1-Cg]'][(tuple(to_explain_statement[0]), tuple(to_explain_statement[1]))]

    prefix = "" if step == 0 else "->"
    alt_x, alt_y = to_explain_statement[0], to_explain_statement[1]
    if all(alt_x - alt_y == 0):
        BDB['[1-Cg]'][(tuple(to_explain_statement[0]), tuple(to_explain_statement[1]))] = True, ("id", None, None, None)
        return True, ("id", None, None, None)

    List_neutral = [i for i in range(len(to_explain_statement[0])) if
                    to_explain_statement[0][i] == to_explain_statement[1][i] == 0 or to_explain_statement[0][i] ==
                    to_explain_statement[1][i] == 1]
    patternRef = [to_explain_statement[0][i] for i in List_neutral]
    if len(List_neutral) == 0:
        BDB['[1-Cg]'][(tuple(to_explain_statement[0]), tuple(to_explain_statement[1]))] = False, (None, None, None, None)

        return False, (
            None, None, None,
            None)  # en réalité si pas de neutre alors CovD0_robust ou Tr_robust ou Tre_robust. Or un 1-Cg verifie d'abord que ceux-ci retournent False

    L_0_1 = [(0, 1)] * len(List_neutral)
    PatternList = list(itert.product(*L_0_1))
    PatternList.remove(tuple(patternRef))
    for pattern in [tuple(patternRef)] + PatternList:
        alt_xr, alt_yr = np.copy(to_explain_statement[0]), np.copy(to_explain_statement[1])
        for i_ in range(len(pattern)):
            alt_xr[List_neutral[i_]] = pattern[i_]
            alt_yr[List_neutral[i_]] = pattern[i_]

        is_explainable_tr_trE_covD0, details = __tr_trE_covD0_robust((alt_xr, alt_yr), PI_statements_List,
                                                                     PI_cofactors_List, step=step, BDB=BDB)
        if is_explainable_tr_trE_covD0:
            _, _, xpl, longueur = details
            BDB['[1-Cg]'][(tuple(to_explain_statement[0]), tuple(to_explain_statement[1]))] = True, (
                prefix + "[1-Cg]", pattern, [chr(ord('a') + neu) for neu in List_neutral], details,
                longueur + (1 if list(pattern) != patternRef else 0))
            return True, (
                prefix + "[1-Cg]", pattern, [chr(ord('a') + neu) for neu in List_neutral], details,
                longueur + (1 if list(pattern) != patternRef else 0))

    BDB['[1-Cg]'][(tuple(to_explain_statement[0]), tuple(to_explain_statement[1]))] = False, (None, None, None, None)
    return False, (None, None, None, None)


def un_Cg_D_robust(to_explain_statement, PI_statements_List, PI_cofactors_List, step=0, BDB=None):
    if BDB is None:
        BDB = dict()

    if (tuple(to_explain_statement[0]), tuple(to_explain_statement[1])) in BDB['[1-CgD]']:
        return BDB['[1-CgD]'][(tuple(to_explain_statement[0]), tuple(to_explain_statement[1]))]

    prefix = "" if step == 0 else "->"
    # identifier les statements du milieu
    identifiant_statement_milieu_List = list()
    for st_selected_ in range(len(PI_statements_List)):
        cofactor_x_z1 = to_explain_statement[0] - PI_statements_List[st_selected_][0]
        cofactor_z2_y = PI_statements_List[st_selected_][1] - to_explain_statement[1]
        if np_computation(cofactor_x_z1, PI_cofactors_List) and np_computation(cofactor_z2_y, PI_cofactors_List):
            identifiant_statement_milieu_List.append(st_selected_)

    if len(identifiant_statement_milieu_List) == 0:
        BDB['[1-CgD]'][(tuple(to_explain_statement[0]), tuple(to_explain_statement[1]))] = False, (None, None, None, None)
        return False, (None, None, None, None)

    Graphe_distances = dict()  # id de x = -1 ; id de y = -2
    Explanations_dict = dict()
    for st_du_milieu in identifiant_statement_milieu_List:
        statement_courant = PI_statements_List[st_du_milieu]
        statement_x_z1_courant = (to_explain_statement[0], statement_courant[0])
        is_cov_explainable, details_cov = cov0_robust(statement_x_z1_courant[0] - statement_x_z1_courant[1],
                                                      PI_cofactors_List, step=step, BDB=BDB)
        if is_cov_explainable:
            Graphe_distances[(-1, st_du_milieu)] = details_cov[-1]
            Explanations_dict[(-1, st_du_milieu)] = details_cov
        else:
            is_un_cg_explainable, details_ = un_Cg_robust(statement_x_z1_courant, PI_statements_List, PI_cofactors_List, step=step, BDB=BDB)

            if is_un_cg_explainable:
                Graphe_distances[(-1, st_du_milieu)] = details_[-1]
                Explanations_dict[(-1, st_du_milieu)] = details_

        statement_z2_y_courant = (statement_courant[1], to_explain_statement[1])
        is_cov_explainable, details_cov = cov0_robust(statement_z2_y_courant[0] - statement_z2_y_courant[1],
                                                      PI_cofactors_List, step=step, BDB=BDB)
        if is_cov_explainable:
            Graphe_distances[(st_du_milieu, -2)] = details_cov[-1]
            Explanations_dict[(st_du_milieu, -2)] = details_cov
        else:
            is_un_cg_explainable, details_ = un_Cg_robust(
                statement_z2_y_courant, PI_statements_List, PI_cofactors_List, step=step, BDB=BDB)
            if is_un_cg_explainable:
                Graphe_distances[(st_du_milieu, -2)] = details_[-1]
                Explanations_dict[(st_du_milieu, -2)] = details_

    # Completer graphe des distances avec les paires de statements du milieu
    for st_du_milieu1 in identifiant_statement_milieu_List:
        for st_du_milieu3 in identifiant_statement_milieu_List:
            if st_du_milieu1 < st_du_milieu3:
                statement1_2 = PI_statements_List[st_du_milieu1]
                statement3_4 = PI_statements_List[st_du_milieu3]

                is_cov_explainable_2_3, details_cov_2_3 = cov0_robust(statement1_2[1] - statement3_4[0],
                                                                      PI_cofactors_List, step=step, BDB=BDB)
                if is_cov_explainable_2_3:
                    Graphe_distances[(st_du_milieu1, st_du_milieu3)] = details_cov_2_3[-1]
                    Explanations_dict[(st_du_milieu1, st_du_milieu3)] = details_cov_2_3
                else:
                    is_rtr_explainable_2_3, details_2_3 = un_Cg_robust(
                        (statement1_2[1], statement3_4[0]), PI_statements_List, PI_cofactors_List, step=step, BDB=BDB)

                    if is_rtr_explainable_2_3:
                        Graphe_distances[(st_du_milieu1, st_du_milieu3)] = details_2_3[-1]
                        Explanations_dict[(st_du_milieu1, st_du_milieu3)] = details_2_3

                is_cov_explainable_4_1, details_cov_4_1 = cov0_robust(statement3_4[1] - statement1_2[0],
                                                                      PI_cofactors_List, step=step, BDB=BDB)
                if is_cov_explainable_4_1:
                    Graphe_distances[(st_du_milieu3, st_du_milieu1)] = details_cov_4_1[-1]
                    Explanations_dict[(st_du_milieu3, st_du_milieu1)] = details_cov_4_1
                else:
                    is_rtr_explainable_4_1, details_4_1 = un_Cg_robust(
                        (statement3_4[1], statement1_2[0]), PI_statements_List, PI_cofactors_List, step=step, BDB=BDB)
                    if is_rtr_explainable_4_1:
                        Graphe_distances[(st_du_milieu3, st_du_milieu1)] = details_4_1[-1]
                        Explanations_dict[(st_du_milieu3, st_du_milieu1)] = details_4_1

    # Calcul du PCC
    model_gurobi = Model("Calcul PCC for explanations")
    model_gurobi.setParam('OutputFlag', False)
    Arcs_Liaisons_Var_Dict = {arc: model_gurobi.addVar(vtype=GRB.BINARY, name=f'b_{arc}') for arc in Graphe_distances}

    # contrainte : exactement un arc sortant du sommet -1
    model_gurobi.addConstr(quicksum([var_b for arc, var_b in Arcs_Liaisons_Var_Dict.items() if arc[0] == -1]) == 1)

    # contraintes : continuité du chemin
    for st_du_milieu in identifiant_statement_milieu_List:
        model_gurobi.addConstr(
            quicksum([var_b for arc, var_b in Arcs_Liaisons_Var_Dict.items() if arc[0] == st_du_milieu]) == quicksum(
                [var_b for arc, var_b in Arcs_Liaisons_Var_Dict.items() if arc[1] == st_du_milieu]))

    model_gurobi.update()
    model_gurobi.setObjective(
        quicksum([(1 + Graphe_distances[arc]) * var_b for arc, var_b in Arcs_Liaisons_Var_Dict.items()]), GRB.MINIMIZE)
    model_gurobi.optimize()

    if model_gurobi.status == GRB.OPTIMAL:
        List_Arcs_PCC = [-1]
        courant = List_Arcs_PCC[-1]
        while courant != -2:
            suiv_courant = -2
            for st_du_milieu in identifiant_statement_milieu_List:
                if (courant, st_du_milieu) in Arcs_Liaisons_Var_Dict and round(
                        Arcs_Liaisons_Var_Dict[(courant, st_du_milieu)].x) == 1:
                    suiv_courant = st_du_milieu
                    List_Arcs_PCC.append(suiv_courant)
            courant = suiv_courant
        List_Arcs_PCC.append(-2)

        explList = list()
        for j in range(1, len(List_Arcs_PCC)):
            explList.append(Explanations_dict[(List_Arcs_PCC[j - 1], List_Arcs_PCC[j])])
            if List_Arcs_PCC[j] != -2:
                explList.append(
                    ("".join([chr(ord('a') + i) for i in range(len(PI_statements_List[List_Arcs_PCC[j]][0])) if
                              PI_statements_List[List_Arcs_PCC[j]][0][i] == 1]),
                     "".join([chr(ord('a') + i) for i in range(len(PI_statements_List[List_Arcs_PCC[j]][1])) if
                              PI_statements_List[List_Arcs_PCC[j]][1][i] == 1])))
            else:
                explList.append(List_Arcs_PCC[j])
        assert explList[-1] == -2

        BDB['[1-CgD]'][(tuple(to_explain_statement[0]), tuple(to_explain_statement[1]))] = True, (prefix + "[1-CgD]", None, explList[:-1], model_gurobi.objVal - 1)
        return True, (prefix + "[1-CgD]", None, explList[:-1], model_gurobi.objVal - 1)

    BDB['[1-CgD]'][(tuple(to_explain_statement[0]), tuple(to_explain_statement[1]))] = False, (None, None, None, None)
    return False, (None, None, None, None)


def deux_Cg_robust(to_explain_statement, PI_statements_List, PI_cofactors_List, step=0, BDB=None):
    if BDB is None:
        BDB = dict()

    if (tuple(to_explain_statement[0]), tuple(to_explain_statement[1])) in BDB['[2-Cg]']:
        return BDB['[2-Cg]'][(tuple(to_explain_statement[0]), tuple(to_explain_statement[1]))]

    prefix = "" if step == 0 else "->"
    alt_x, alt_y = to_explain_statement[0], to_explain_statement[1]
    if all(alt_x - alt_y == 0):
        BDB['[2-Cg]'][(tuple(to_explain_statement[0]), tuple(to_explain_statement[1]))] = True, ("id", None, None, None)
        return True, ("id", None, None, None)

    List_neutral = [i for i in range(len(to_explain_statement[0])) if
                    to_explain_statement[0][i] == to_explain_statement[1][i] == 0 or to_explain_statement[0][i] ==
                    to_explain_statement[1][i] == 1]
    patternRef = [to_explain_statement[0][i] for i in List_neutral]
    if len(List_neutral) == 0:
        BDB['[2-Cg]'][(tuple(to_explain_statement[0]), tuple(to_explain_statement[1]))] = False, (None, None, None, None)
        return False, (
            None, None, None,
            None)  # raisonnement identique à un_Cg_robuste

    L_0_1 = [(0, 1)] * len(List_neutral)
    PatternList = list(itert.product(*L_0_1))
    PatternList.remove(tuple(patternRef))
    for pattern in [tuple(patternRef)] + PatternList:
        alt_xr, alt_yr = np.copy(to_explain_statement[0]), np.copy(to_explain_statement[1])
        for i_ in range(len(pattern)):
            alt_xr[List_neutral[i_]] = pattern[i_]
            alt_yr[List_neutral[i_]] = pattern[i_]

        is_explainable_one_cg_d, details = un_Cg_D_robust((alt_xr, alt_yr), PI_statements_List,
                                                          PI_cofactors_List, step=step, BDB=BDB)
        if is_explainable_one_cg_d:
            _, _, xpl, longueur = details
            BDB['[2-Cg]'][(tuple(to_explain_statement[0]), tuple(to_explain_statement[1]))] = True, (
                prefix + "[2-Cg]", pattern, [chr(ord('a') + neu) for neu in List_neutral], details,
                longueur + (1 if list(pattern) != patternRef else 0))
            return True, (
                prefix + "[2-Cg]", pattern, [chr(ord('a') + neu) for neu in List_neutral], details,
                longueur + (1 if list(pattern) != patternRef else 0))

    BDB['[2-Cg]'][(tuple(to_explain_statement[0]), tuple(to_explain_statement[1]))] = False, (None, None, None, None)
    return False, (None, None, None, None)


def __cov_un_cg_deux_cg(to_explain_statement, PI_statements_List, PI_cofactors_List, step=0, BDB=None):
    if BDB is None:
        BDB = dict()

    is_cov, details_cov = cov0_robust(to_explain_statement[0] - to_explain_statement[1], PI_cofactors_List, step=step, BDB=BDB)
    if is_cov:
        return True, details_cov
    is_un_cg, details_un_cg = un_Cg_robust(to_explain_statement, PI_statements_List, PI_cofactors_List, step=step, BDB=BDB)
    if is_un_cg:
        return True, details_un_cg
    return deux_Cg_robust(to_explain_statement, PI_statements_List, PI_cofactors_List, step=step, BDB=BDB)


def deux_Cg_D_robust(to_explain_statement, PI_statements_List, PI_cofactors_List, step=0, BDB=None):
    if BDB is None:
        BDB = dict()

    if (tuple(to_explain_statement[0]), tuple(to_explain_statement[1])) in BDB['[2-CgD]']:
        return BDB['[2-CgD]'][(tuple(to_explain_statement[0]), tuple(to_explain_statement[1]))]

    prefix = "" if step == 0 else "->"
    # identifier les statements du milieu
    identifiant_statement_milieu_List = list()
    for st_selected_ in range(len(PI_statements_List)):
        cofactor_x_z1 = to_explain_statement[0] - PI_statements_List[st_selected_][0]
        cofactor_z2_y = PI_statements_List[st_selected_][1] - to_explain_statement[1]
        if np_computation(cofactor_x_z1, PI_cofactors_List) and np_computation(cofactor_z2_y, PI_cofactors_List):
            identifiant_statement_milieu_List.append(st_selected_)

    if len(identifiant_statement_milieu_List) == 0:
        BDB['[2-CgD]'][(tuple(to_explain_statement[0]), tuple(to_explain_statement[1]))] = False, (None, None, None, None)
        return False, (None, None, None, None)

    Graphe_distances = dict()  # id de x = -1 ; id de y = -2
    Explanations_dict = dict()
    for st_du_milieu in identifiant_statement_milieu_List:
        statement_courant = PI_statements_List[st_du_milieu]
        statement_x_z1_courant = (to_explain_statement[0], statement_courant[0])
        is_rtr_explainable_, details_ = __cov_un_cg_deux_cg(statement_x_z1_courant, PI_statements_List,
                                                            PI_cofactors_List, step=step, BDB=BDB)

        if is_rtr_explainable_:
            Graphe_distances[(-1, st_du_milieu)] = details_[-1]
            Explanations_dict[(-1, st_du_milieu)] = details_

        statement_z2_y_courant = (statement_courant[1], to_explain_statement[1])
        is_rtr_explainable_, details_ = __cov_un_cg_deux_cg(
            statement_z2_y_courant, PI_statements_List, PI_cofactors_List, step=step, BDB=BDB)
        if is_rtr_explainable_:
            Graphe_distances[(st_du_milieu, -2)] = details_[-1]
            Explanations_dict[(st_du_milieu, -2)] = details_

    # Completer graphe des distances avec les paires de statements du milieu
    for st_du_milieu1 in identifiant_statement_milieu_List:
        for st_du_milieu3 in identifiant_statement_milieu_List:
            if st_du_milieu1 < st_du_milieu3:
                statement1_2 = PI_statements_List[st_du_milieu1]
                statement3_4 = PI_statements_List[st_du_milieu3]

                is_rtr_explainable_2_3, details_2_3 = __cov_un_cg_deux_cg(
                    (statement1_2[1], statement3_4[0]), PI_statements_List, PI_cofactors_List, step=step, BDB=BDB)

                if is_rtr_explainable_2_3:
                    Graphe_distances[(st_du_milieu1, st_du_milieu3)] = details_2_3[-1]
                    Explanations_dict[(st_du_milieu1, st_du_milieu3)] = details_2_3

                is_rtr_explainable_4_1, details_4_1 = __cov_un_cg_deux_cg(
                    (statement3_4[1], statement1_2[0]), PI_statements_List, PI_cofactors_List, step=step, BDB=BDB)
                if is_rtr_explainable_4_1:
                    Graphe_distances[(st_du_milieu3, st_du_milieu1)] = details_4_1[-1]
                    Explanations_dict[(st_du_milieu3, st_du_milieu1)] = details_4_1

    # Calcul du PCC
    model_gurobi = Model("Calcul PCC for explanations")
    model_gurobi.setParam('OutputFlag', False)
    Arcs_Liaisons_Var_Dict = {arc: model_gurobi.addVar(vtype=GRB.BINARY, name=f'b_{arc}') for arc in Graphe_distances}

    # contrainte : exactement un arc sortant du sommet -1
    model_gurobi.addConstr(quicksum([var_b for arc, var_b in Arcs_Liaisons_Var_Dict.items() if arc[0] == -1]) == 1)

    # contraintes : continuité du chemin
    for st_du_milieu in identifiant_statement_milieu_List:
        model_gurobi.addConstr(
            quicksum([var_b for arc, var_b in Arcs_Liaisons_Var_Dict.items() if arc[0] == st_du_milieu]) == quicksum(
                [var_b for arc, var_b in Arcs_Liaisons_Var_Dict.items() if arc[1] == st_du_milieu]))

    model_gurobi.update()
    model_gurobi.setObjective(
        quicksum([(1 + Graphe_distances[arc]) * var_b for arc, var_b in Arcs_Liaisons_Var_Dict.items()]), GRB.MINIMIZE)
    model_gurobi.optimize()

    if model_gurobi.status == GRB.OPTIMAL:
        List_Arcs_PCC = [-1]
        courant = List_Arcs_PCC[-1]
        while courant != -2:
            suiv_courant = -2
            for st_du_milieu in identifiant_statement_milieu_List:
                if (courant, st_du_milieu) in Arcs_Liaisons_Var_Dict and round(
                        Arcs_Liaisons_Var_Dict[(courant, st_du_milieu)].x) == 1:
                    suiv_courant = st_du_milieu
                    List_Arcs_PCC.append(suiv_courant)
            courant = suiv_courant
        List_Arcs_PCC.append(-2)

        explList = list()
        for j in range(1, len(List_Arcs_PCC)):
            explList.append(Explanations_dict[(List_Arcs_PCC[j - 1], List_Arcs_PCC[j])])
            if List_Arcs_PCC[j] != -2:
                explList.append(
                    ("".join([chr(ord('a') + i) for i in range(len(PI_statements_List[List_Arcs_PCC[j]][0])) if
                              PI_statements_List[List_Arcs_PCC[j]][0][i] == 1]),
                     "".join([chr(ord('a') + i) for i in range(len(PI_statements_List[List_Arcs_PCC[j]][1])) if
                              PI_statements_List[List_Arcs_PCC[j]][1][i] == 1])))
            else:
                explList.append(List_Arcs_PCC[j])
        assert explList[-1] == -2

        BDB['[2-CgD]'][(tuple(to_explain_statement[0]), tuple(to_explain_statement[1]))] = True, (prefix + "[2-CgD]", None, explList[:-1], model_gurobi.objVal - 1)
        return True, (prefix + "[2-CgD]", None, explList[:-1], model_gurobi.objVal - 1)

    BDB['[2-CgD]'][(tuple(to_explain_statement[0]), tuple(to_explain_statement[1]))] = False, (None, None, None, None)
    return False, (None, None, None, None)
