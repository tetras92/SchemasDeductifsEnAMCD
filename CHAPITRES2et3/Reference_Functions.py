import numpy as np
from gurobipy import *

import itertools as itert
from Utils import check_int_List_1, check_int_List_2, integer_to_bin

from CBTOprocessing import flat_CBTO_formated_for_OfflineSimulator, Tn


def _doPa(covector_to_check, step=0, BDB=None):
    if BDB is None:
        BDB = dict()

    prefix = "" if step == 0 else "->"
    if all(covector_to_check >= 0) and not (all(covector_to_check == 0)):
        return True, (prefix + "[doPa]", None, None, 1)

    return False, (None, None, None, None)


def tr(to_explain_statement, PI_statements_List, W, step=0, BDB=None):
    if BDB is None:
        BDB = dict()
    prefix = "" if step == 0 else "->"

    alt_x, alt_y = to_explain_statement[0], to_explain_statement[1]
    if (tuple(alt_x), tuple(alt_y)) in BDB['[tr]']:
        return BDB['[tr]'][(tuple(alt_x), tuple(alt_y))]

    def pareto_domine(dominant_array, dominated_array):
        return all(dominant_array - dominated_array >= 0) and not (
            all(dominant_array - dominated_array == 0))  #

    Successeur_Dict = {st_k: [st_k_prime for st_k_prime in range(len(PI_statements_List)) if
                              np.all(PI_statements_List[st_k][1] == PI_statements_List[st_k_prime][0])
                              or pareto_domine(PI_statements_List[st_k][1], PI_statements_List[st_k_prime][0])] + (
                                 [-2] if pareto_domine(PI_statements_List[st_k][1], alt_y) or np.all(
                                     PI_statements_List[st_k][1] == alt_y) else list()) for st_k in
                       range(len(PI_statements_List))}
    Successeur_Dict[-1] = [st_k for st_k in range(len(PI_statements_List)) if
                           pareto_domine(alt_x, PI_statements_List[st_k][0]) or np.all(
                               alt_x == PI_statements_List[st_k][0])]

    Successeur_Dict[-2] = list()

    model_gurobi = Model("Calcul PCC for transitive scheme explanations")
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
        # print(PresenceDominance)
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
        details = (prefix + "[tr]",
                   Chaine_des_st_k[1:-1],
                   Chaine_a_retourner,
                   len(Chaine_des_st_k[1:-1]))
        BDB['[tr]'][(tuple(alt_x), tuple(alt_y))] = True, details
        return True, details

    BDB['[tr]'][(tuple(alt_x), tuple(alt_y))] = (False, (None, None, None, None))
    return False, (None, None, None, None)


def cg(to_explain_statement, PI_statements_List, W, step=0, BDB=None):
    PI_covectors_List = [x_arr - y_arr for x_arr, y_arr in PI_statements_List]
    if BDB is None:
        BDB = dict()
    alt_x, alt_y = to_explain_statement[0], to_explain_statement[1]
    if (tuple(alt_x), tuple(alt_y)) in BDB['[cg]']:
        return BDB['[cg]'][(tuple(alt_x), tuple(alt_y))]

    if all(alt_x - alt_y == 0):
        return True, ("id", None, None, 0)

    is_dopa, details_dopa = _doPa(alt_x - alt_y, step, BDB)  # à jour ce 06/10
    if is_dopa:
        BDB['[cg]'][(tuple(alt_x), tuple(alt_y))] = True, details_dopa
        return True, details_dopa

    prefix = "" if step == 0 else "->"
    covector_to_check = alt_x - alt_y
    for k in range(len(PI_covectors_List)):
        cov_pi = PI_covectors_List[k]
        if all(covector_to_check - cov_pi == 0):
            details = (prefix + "[cg]", ("".join(
                [chr(ord('a') + i) for i in range(len(PI_statements_List[k][0])) if PI_statements_List[k][0][i] == 1]),
                                         "".join(
                                             [chr(ord('a') + i) for i in range(len(PI_statements_List[k][1])) if
                                              PI_statements_List[k][1][i] == 1])), None, 1)
            BDB['[cg]'][(tuple(alt_x), tuple(alt_y))] = True, details

            return True, details

    BDB['[cg]'][(tuple(alt_x), tuple(alt_y))] = False, (None, None, None, None)
    return False, (None, None, None, None)


def cov(to_explain_statement, PI_statements_List, W, step=0, BDB=None):
    if BDB is None:
        BDB = dict()
    prefix = "" if step == 0 else "->"
    alt_x, alt_y = to_explain_statement[0], to_explain_statement[1]
    if (tuple(alt_x), tuple(alt_y)) in BDB['[cov]']:
        return BDB['[cov]'][(tuple(alt_x), tuple(alt_y))]

    alt_x, alt_y = to_explain_statement[0], to_explain_statement[1]

    SigmaP = {i_ for i_ in range(len(alt_x)) if alt_x[i_] == 1}
    SigmaD = {i_ for i_ in range(len(alt_y)) if alt_y[i_] == 1}

    model_gurobi = Model("Calcul Couverture general")
    model_gurobi.setParam('OutputFlag', False)

    BVar_Dict = {j: model_gurobi.addVar(vtype=GRB.BINARY, name=f'b_{j}') for
                 j in range(len(PI_statements_List))}

    # C1 : aucune aide pour les Pros
    for i in set(range(len(alt_x))) - SigmaP:
        model_gurobi.addConstr(
            quicksum([BVar_Dict[j] for j in range(len(PI_statements_List)) if PI_statements_List[j][0][i] == 1]) == 0)

    # C2
    for i in SigmaP:
        model_gurobi.addConstr(
            quicksum([BVar_Dict[j] for j in range(len(PI_statements_List)) if PI_statements_List[j][0][i] == 1]) <= 1)

    # C3
    for i in SigmaD:
        model_gurobi.addConstr(
            quicksum([BVar_Dict[j] for j in range(len(PI_statements_List)) if PI_statements_List[j][1][i] == 1]) == 1)

    # C4
    for i in set(range(len(alt_y))) - SigmaD:  #
        model_gurobi.addConstr(
            quicksum([BVar_Dict[j] for j in range(len(PI_statements_List)) if PI_statements_List[j][1][i] == 1]) <= 1)

    # C5
    model_gurobi.addConstr(quicksum(BVar_Dict.values()) >= 1)  # non inDispensable

    model_gurobi.update()

    model_gurobi.setObjective(quicksum(BVar_Dict.values()), GRB.MINIMIZE)
    model_gurobi.optimize()
    if model_gurobi.status == GRB.OPTIMAL:
        Xprime = {i for j in range(len(PI_statements_List)) if BVar_Dict[j].x == 1 for i in
                  range(len(alt_x)) if PI_statements_List[j][0][i] == 1}

        Yprime = {i for j in range(len(PI_statements_List)) if BVar_Dict[j].x == 1 for i in
                  range(len(alt_x)) if PI_statements_List[j][1][i] == 1}

        is_SigmaP_Subset_Xprime = len(SigmaP - Xprime) > 0

        is_Yprime_Subset_SigmaD = len(Yprime - SigmaD) > 0

        L = list()
        for j in range(len(PI_statements_List)):
            if BVar_Dict[j].x == 1:
                set_uj = {i for i in range(len(alt_x)) if PI_statements_List[j][0][i] == 1}
                set_vj = {i for i in range(len(alt_x)) if
                          PI_statements_List[j][1][i] == 1}
                L.append(("".join([chr(ord('a') + i) for i in set_uj]), "".join([chr(ord('a') + i) for i in set_vj])))

        if is_SigmaP_Subset_Xprime:
            L = ['D'] + L

        if is_Yprime_Subset_SigmaD:
            L += ['D']

        details = (prefix + "[cov]", None, L, len(L))
        BDB['[cov]'][(tuple(alt_x), tuple(alt_y))] = (True, details)
        return True, details

    BDB['[cov]'][(tuple(alt_x), tuple(alt_y))] = False, (None, None, None, None)
    return False, (None, None, None, None)


def cg_generique(to_explain_statement, PI_statements_List, W, corresponding_function, step=0, BDB=None):
    if BDB is None:
        BDB = dict()

    name = None
    if corresponding_function is tr:
        name = '[cg/tr]'
    elif corresponding_function is cov:
        name = '[cg/cov]'
    elif corresponding_function is tr_cov:
        name = '[cg/tr/cov]'
    elif corresponding_function is cov_tr:
        name = '[cg/cov/tr]'
    else:
        raise Exception("cg - generique")

    if (tuple(to_explain_statement[0]), tuple(to_explain_statement[1])) in BDB[name]:
        return BDB[name][(tuple(to_explain_statement[0]), tuple(to_explain_statement[1]))]

    prefix = "" if step == 0 else "->"

    List_neutral = [i for i in range(len(to_explain_statement[0])) if
                    to_explain_statement[0][i] == to_explain_statement[1][i] == 0 or to_explain_statement[0][i] ==
                    to_explain_statement[1][i] == 1]
    patternRef = [to_explain_statement[0][i] for i in List_neutral]
    if len(List_neutral) == 0:
        BDB[name][(tuple(to_explain_statement[0]), tuple(to_explain_statement[1]))] = False, (
            None, None, None, None)
        return False, (
            None, None, None,
            None)  #

    L_0_1 = [(0, 1)] * len(List_neutral)
    PatternList = list(itert.product(*L_0_1))
    PatternList.remove(tuple(patternRef))
    for pattern in [tuple(patternRef)] + PatternList:
        alt_xr, alt_yr = np.copy(to_explain_statement[0]), np.copy(to_explain_statement[1])
        for i_ in range(len(pattern)):
            alt_xr[List_neutral[i_]] = pattern[i_]
            alt_yr[List_neutral[i_]] = pattern[i_]

        is_intermediate_explainable, details = corresponding_function((alt_xr, alt_yr), PI_statements_List, W,
                                                                      step=step, BDB=BDB)
        if is_intermediate_explainable:
            _, _, xpl, longueur = details
            BDB[name][(tuple(to_explain_statement[0]), tuple(to_explain_statement[1]))] = True, (
                prefix + name, pattern, [chr(ord('a') + neu) for neu in List_neutral], details,
                longueur + (1 if list(pattern) != patternRef else 0))
            return True, (
                prefix + name, pattern, [chr(ord('a') + neu) for neu in List_neutral], details,
                longueur + (1 if list(pattern) != patternRef else 0))

    BDB[name][(tuple(to_explain_statement[0]), tuple(to_explain_statement[1]))] = False, (None, None, None, None)
    return False, (None, None, None, None)


def cg_cov(to_explain_statement, PI_statements_List, W, step=0, BDB=None):
    alt_x, alt_y = to_explain_statement[0], to_explain_statement[1]
    if all(alt_x - alt_y == 0):
        return True, ("id", None, None, 0)
    return cg_generique(to_explain_statement, PI_statements_List, W, cov, step, BDB)


def cg_tr(to_explain_statement, PI_statements_List, W, step=0, BDB=None):
    return cg_generique(to_explain_statement, PI_statements_List, W, tr, step, BDB)


def tr_generique(to_explain_statement, PI_statements_List, W, corresponding_function, step=0, BDB=None):
    if BDB is None:
        BDB = dict()

    name = None
    if corresponding_function is cg:
        name = '[tr/cg]'
    elif corresponding_function is cov:
        name = '[tr/cov]'
    elif corresponding_function is cg_cov:
        name = '[tr/cg/cov]'
    elif corresponding_function is cov_cg:
        name = '[tr/cov/cg]'
    else:
        raise Exception("tr - generique", corresponding_function)

    if (tuple(to_explain_statement[0]), tuple(to_explain_statement[1])) in BDB[name]:
        return BDB[name][(tuple(to_explain_statement[0]), tuple(to_explain_statement[1]))]

    prefix = "" if step == 0 else "->"
    # identifier les statements du milieu
    identifiant_statement_milieu_List = list()
    X_score = sum([to_explain_statement[0][i] * w_i for i, w_i in W.items()])
    Y_score = sum([to_explain_statement[1][i] * w_i for i, w_i in W.items()])

    for st_selected_ in range(len(PI_statements_List)):
        z1_score = sum([PI_statements_List[st_selected_][0][i] * w_i for i, w_i in W.items()])
        z2_score = sum([PI_statements_List[st_selected_][1][i] * w_i for i, w_i in W.items()])
        if X_score >= z1_score >= z2_score >= Y_score:
            identifiant_statement_milieu_List.append(st_selected_)

    Graphe_distances = dict()  # id de x = -1 ; id de y = -2
    Explanations_dict = dict()

    is_direct_explainable, details_direct = corresponding_function((to_explain_statement[0], to_explain_statement[1]),
                                                                   PI_statements_List, W, step, BDB)
    if is_direct_explainable:
        Graphe_distances[(-1, -2)] = details_direct[-1]
        Explanations_dict[(-1, -2)] = details_direct

    for st_du_milieu in identifiant_statement_milieu_List:
        statement_courant = PI_statements_List[st_du_milieu]
        statement_x_z1_courant = (to_explain_statement[0], statement_courant[0])
        is_intermediate_explainable, details_cov = corresponding_function(
            (statement_x_z1_courant[0], statement_x_z1_courant[1]),
            PI_statements_List, W,
            step=step, BDB=BDB)
        if is_intermediate_explainable:
            Graphe_distances[(-1, st_du_milieu)] = details_cov[-1]
            Explanations_dict[(-1, st_du_milieu)] = details_cov

        statement_z2_y_courant = (statement_courant[1], to_explain_statement[1])
        is_intermediate_explainable, details_cov = corresponding_function(
            (statement_z2_y_courant[0], statement_z2_y_courant[1]),
            PI_statements_List, W,
            step=step, BDB=BDB)
        if is_intermediate_explainable:
            Graphe_distances[(st_du_milieu, -2)] = details_cov[-1]
            Explanations_dict[(st_du_milieu, -2)] = details_cov

    # Completer graphe des distances avec les paires de statements du milieu
    for st_du_milieu1 in identifiant_statement_milieu_List:
        for st_du_milieu3 in identifiant_statement_milieu_List:
            if st_du_milieu1 < st_du_milieu3:
                statement1_2 = PI_statements_List[st_du_milieu1]
                statement3_4 = PI_statements_List[st_du_milieu3]

                is_intermediate_explainable_2_3, details_cov_2_3 = corresponding_function(
                    (statement1_2[1], statement3_4[0]),
                    PI_statements_List, W,
                    step=step, BDB=BDB)
                if is_intermediate_explainable_2_3:
                    Graphe_distances[(st_du_milieu1, st_du_milieu3)] = details_cov_2_3[-1]
                    Explanations_dict[(st_du_milieu1, st_du_milieu3)] = details_cov_2_3

                is_intermediate_explainable_4_1, details_cov_4_1 = corresponding_function(
                    (statement3_4[1], statement1_2[0]),
                    PI_statements_List, W,
                    step=step, BDB=BDB)
                if is_intermediate_explainable_4_1:
                    Graphe_distances[(st_du_milieu3, st_du_milieu1)] = details_cov_4_1[-1]
                    Explanations_dict[(st_du_milieu3, st_du_milieu1)] = details_cov_4_1

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

        BDB[name][(tuple(to_explain_statement[0]), tuple(to_explain_statement[1]))] = True, (
            prefix + name, None, explList[:-1], model_gurobi.objVal - 1)
        return True, (prefix + name, None, explList[:-1], model_gurobi.objVal - 1)

    BDB[name][(tuple(to_explain_statement[0]), tuple(to_explain_statement[1]))] = False, (None, None, None, None)
    return False, (None, None, None, None)


def tr_cov(to_explain_statement, PI_statements_List, W, step=0, BDB=None):
    return tr_generique(to_explain_statement, PI_statements_List, W, cov, step, BDB)


def tr_cg(to_explain_statement, PI_statements_List, W, step=0, BDB=None):
    return tr_generique(to_explain_statement, PI_statements_List, W, cg, step, BDB)


def _cov_incomplet(to_explain_statement, PI_statements_List, W, step=0, BDB=None):
    if BDB is None:
        BDB = dict()
    prefix = "" if step == 0 else "->"
    alt_x, alt_y = to_explain_statement[0], to_explain_statement[1]

    SigmaP = {i_ for i_ in range(len(alt_x)) if alt_x[i_] == 1}
    SigmaD = {i_ for i_ in range(len(alt_y)) if alt_y[i_] == 1}

    model_gurobi = Model("Calcul Couverture general")
    model_gurobi.setParam('OutputFlag', False)

    BVar_Dict = {j: model_gurobi.addVar(vtype=GRB.BINARY, name=f'b_{j}') for
                 j in range(len(PI_statements_List))}

    X_item_Var = {i_: model_gurobi.addVar(vtype=GRB.BINARY, name=f'xv_{i_}') for i_ in SigmaP}
    Y_item_Var = {i_: model_gurobi.addVar(vtype=GRB.BINARY, name=f'yv_{i_}') for i_ in SigmaD}

    # C1 : aucune aide pour les Pros
    for i in set(range(len(alt_x))) - SigmaP:
        model_gurobi.addConstr(
            quicksum([BVar_Dict[j] for j in range(len(PI_statements_List)) if PI_statements_List[j][0][i] == 1]) == 0)

    # C2
    for i in SigmaP:
        model_gurobi.addConstr(
            quicksum([BVar_Dict[j] for j in range(len(PI_statements_List)) if PI_statements_List[j][0][i] == 1]) ==
            X_item_Var[i])

    # C3
    for i in set(range(len(alt_y))) - SigmaD:
        model_gurobi.addConstr(
            quicksum([BVar_Dict[j] for j in range(len(PI_statements_List)) if PI_statements_List[j][1][
                i] == 1]) <= 1)

    # C4
    for i in SigmaD:
        model_gurobi.addConstr(
            quicksum([BVar_Dict[j] for j in range(len(PI_statements_List)) if PI_statements_List[j][1][i] == 1]) ==
            Y_item_Var[i])

    # C5
    model_gurobi.addConstr(quicksum([(1 - X_item_Var[i_]) * W[i_] for i_ in X_item_Var]) >= quicksum(
        [(1 - Y_item_Var[i_]) * W[i_] for i_ in Y_item_Var]) + 1)

    # C6
    model_gurobi.addConstr(quicksum(BVar_Dict.values()) >= 1)  # Corriger aujourd'hui 13/09/23

    model_gurobi.update()

    model_gurobi.setObjective(
        quicksum(X_item_Var.values()) + quicksum(Y_item_Var.values()),
        GRB.MINIMIZE)

    model_gurobi.optimize()
    if model_gurobi.status == GRB.OPTIMAL:

        ElementsDeXRecuperes = {i_ for i_ in SigmaP if X_item_Var[i_].x == 1}
        ElementsDeYRecuperes = {i_ for i_ in SigmaD if Y_item_Var[i_].x == 1}
        # print(ElementsDeXRecuperes, ElementsDeYRecuperes)
        L = list()
        for j in range(len(PI_statements_List)):
            if BVar_Dict[j].x == 1:
                set_uj = {i for i in range(len(alt_x)) if PI_statements_List[j][0][i] == 1}

                set_vj = {i for i in range(len(alt_x)) if
                          PI_statements_List[j][1][i] == 1}
                L.append(("".join([chr(ord('a') + i) for i in set_uj]), "".join([chr(ord('a') + i) for i in set_vj])))

        remaining_pwc = np.zeros(len(alt_x)), np.zeros(len(alt_x))
        for i in SigmaP - ElementsDeXRecuperes:
            remaining_pwc[0][i] = 1

        for i in SigmaD - ElementsDeYRecuperes:
            remaining_pwc[1][i] = 1

        return L + [remaining_pwc]

    return None


def cov_generique(to_explain_statement, PI_statements_List, W, corresponding_function, step=0, BDB=None):
    if BDB is None:
        BDB = dict()

    name = None
    if corresponding_function is cg:
        name = '[cov/cg]'
    elif corresponding_function is tr:
        name = '[cov/tr]'
    elif corresponding_function is cg_tr:
        name = '[cov/cg/tr]'
    elif corresponding_function is tr_cg:
        name = '[cov/tr/cg]'
    else:
        raise Exception("cov - generique")

    if (tuple(to_explain_statement[0]), tuple(to_explain_statement[1])) in BDB[name]:
        return BDB[name][(tuple(to_explain_statement[0]), tuple(to_explain_statement[1]))]

    prefix = "" if step == 0 else "->"
    result_cov_incomplet = _cov_incomplet(to_explain_statement, PI_statements_List, W, step=0, BDB=None)
    if result_cov_incomplet is None:
        BDB[name][(tuple(to_explain_statement[0]), tuple(to_explain_statement[1]))] = False, (None, None, None, None)
        return False, (None, None, None, None)

    is_intermediate_explainable, details = corresponding_function(result_cov_incomplet[-1], PI_statements_List, W,
                                                                  step, BDB)
    if not is_intermediate_explainable:
        BDB[name][(tuple(to_explain_statement[0]), tuple(to_explain_statement[1]))] = False, (None, None, None, None)
        return False, (None, None, None, None)

    details = (prefix + name, result_cov_incomplet[:-1], details, details[-1] + len(result_cov_incomplet[:-1]))
    BDB[name][(tuple(to_explain_statement[0]), tuple(to_explain_statement[1]))] = True, details
    return True, details


def cov_tr(to_explain_statement, PI_statements_List, W, step=0, BDB=None):
    return cov_generique(to_explain_statement, PI_statements_List, W, tr, step, BDB)


def cov_cg(to_explain_statement, PI_statements_List, W, step=0, BDB=None):
    alt_x, alt_y = to_explain_statement[0], to_explain_statement[1]
    if all(alt_x - alt_y == 0):
        return True, ("id", None, None, 0)
    return cov_generique(to_explain_statement, PI_statements_List, W, cg, step, BDB)


###############################################################

def generate_set_A(m, size_A=10):
    setA_as_integers = tuple(*np.random.choice(range(1, 2 ** m - 1), size=(1, size_A), replace=False))
    # while not (check_int_List_1(setA_as_integers, m) and check_int_List_2(setA_as_integers, m)):
    while not (check_int_List_1(setA_as_integers, m)):
        setA_as_integers = tuple(*np.random.choice(range(1, 2 ** m - 1), size=(1, size_A), replace=False))

    setA_as_binary_np_arrays = [np.array([int(binary_digit) for binary_digit in integer_to_bin(int(alt_int_val), m)])
                                for alt_int_val in setA_as_integers]

    return setA_as_binary_np_arrays


def generate_set_AR(m, size_AR=10):
    setAR_as_integers = tuple(*np.random.choice(range(1, 2 ** m - 1), size=(1, size_AR), replace=False))
    while not (check_int_List_1(setAR_as_integers, m) and check_int_List_2(setAR_as_integers, m)):
        # while not (check_int_List_1(setAR_as_integers, m)):
        setAR_as_integers = tuple(*np.random.choice(range(1, 2 ** m - 1), size=(1, size_AR), replace=False))

    setAR_as_binary_np_arrays = [np.array([int(binary_digit) for binary_digit in integer_to_bin(int(alt_int_val), m)])
                                 for alt_int_val in setAR_as_integers]

    return setAR_as_binary_np_arrays


def generate_flat_order(filename, m):
    L = flat_CBTO_formated_for_OfflineSimulator(filename, m)
    L_as_bins = [integer_to_bin(int(alt_int_val), m)[::-1] for alt_int_val in L]
    L_as_sets = [{i for i in range(len(thebin)) if thebin[i] == '1'} for thebin in L_as_bins]
    L_as_str = ["".join([str(elm) for elm in sorted(theset)]) for theset in L_as_sets]

    return L_as_str


def generate_critical_pairs_and_Tm(filename, m):
    L = flat_CBTO_formated_for_OfflineSimulator(filename, m)
    L_as_bins = [integer_to_bin(int(alt_int_val), m)[::-1] for alt_int_val in L]
    L_as_sets = [{i for i in range(len(thebin)) if thebin[i] == '1'} for thebin in L_as_bins]
    TmTransformed = [({p - 1 for p in pair[0]}, {p - 1 for p in pair[1]}) for pair in Tn(m)]
    # print(L_as_sets)
    STn = [(min(pair, key=lambda x: L_as_sets.index(x)),
            max(pair, key=lambda x: L_as_sets.index(x))) for pair in TmTransformed]

    CriticalPair = [(L_as_sets[j], L_as_sets[j + 1]) for j in range(len(L_as_sets) - 1) if
                    (L_as_sets[j], L_as_sets[j + 1]) in STn]  # contigu et disjoints
    return CriticalPair, STn


def generate_critical_pairs_and_Atomic(filename, m):
    L = flat_CBTO_formated_for_OfflineSimulator(filename, m)
    L_as_bins = [integer_to_bin(int(alt_int_val), m)[::-1] for alt_int_val in L]
    L_as_sets = [{i for i in range(len(thebin)) if thebin[i] == '1'} for thebin in L_as_bins]
    TmTransformed = [({p - 1 for p in pair[0]}, {p - 1 for p in pair[1]}) for pair in Tn(m)]
    # print(L_as_sets)
    STn = [(min(pair, key=lambda x: L_as_sets.index(x)),
            max(pair, key=lambda x: L_as_sets.index(x))) for pair in TmTransformed]

    Atomic = [(min(pair, key=lambda x: L_as_sets.index(x)),
               max(pair, key=lambda x: L_as_sets.index(x))) for pair in TmTransformed if
              len(pair[0]) == 1 or len(pair[1]) == 1]

    CriticalPair = [(L_as_sets[j], L_as_sets[j + 1]) for j in range(len(L_as_sets) - 1) if
                    (L_as_sets[j], L_as_sets[j + 1]) in STn]  # contigu et disjoints
    return CriticalPair, Atomic

