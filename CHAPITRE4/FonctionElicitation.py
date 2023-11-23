from Tools import CONSTRAINTSFEASIBILITYTOL
from gurobipy import *
import numpy as np
import itertools as iter_


def is_2Order_GMS_with_PoolSolutions(to_explain_statement, PI_statements_List, PI_cofactors_List, poolSolutions=5):
    Graphe_distances = dict()  # id de x = -1 ; id de y = -2
    Pros_Dict = dict()
    Cons_Dict = dict()

    # pour toute paire de statements
    for k, k_prime in iter_.combinations(range(len(PI_statements_List)), 2):
        st_k_plus, st_k_moins = PI_statements_List[k][0], PI_statements_List[k][1]
        st_k_prime_plus, st_k_prime_moins = PI_statements_List[k_prime][0], PI_statements_List[k_prime][1]

        k_k_prime_cofactor = st_k_moins - st_k_prime_plus
        pro_k_k_prime = {i_ for i_ in range(len(k_k_prime_cofactor)) if k_k_prime_cofactor[i_] == 1}
        con_k_k_prime = {j_ for j_ in range(len(k_k_prime_cofactor)) if k_k_prime_cofactor[j_] == -1}
        Pros_Dict[(k, k_prime)] = pro_k_k_prime
        Cons_Dict[(k, k_prime)] = con_k_k_prime
        Graphe_distances[(k, k_prime)] = 0 if len(con_k_k_prime) > len(pro_k_k_prime) else (
                len(con_k_k_prime) + (1 if len(con_k_k_prime) < len(pro_k_k_prime) else 0))

        k_prime_k_cofactor = st_k_prime_moins - st_k_plus
        pro_k_prime_k = {i_ for i_ in range(len(k_prime_k_cofactor)) if k_prime_k_cofactor[i_] == 1}
        con_k_prime_k = {j_ for j_ in range(len(k_prime_k_cofactor)) if k_prime_k_cofactor[j_] == -1}
        Pros_Dict[(k_prime, k)] = pro_k_prime_k
        Cons_Dict[(k_prime, k)] = con_k_prime_k
        Graphe_distances[(k_prime, k)] = 0 if len(con_k_prime_k) > len(pro_k_prime_k) else (
                len(con_k_prime_k) + (1 if len(con_k_prime_k) < len(pro_k_prime_k) else 0))

    alt_x, alt_y = to_explain_statement
    for k in range(len(PI_statements_List)):
        st = PI_statements_List[k]
        x_st_cofactor = alt_x - st[0]
        pro_x_st = {i_ for i_ in range(len(x_st_cofactor)) if x_st_cofactor[i_] == 1}
        con_x_st = {j_ for j_ in range(len(x_st_cofactor)) if x_st_cofactor[j_] == -1}
        Pros_Dict[(-1, k)] = pro_x_st
        Cons_Dict[(-1, k)] = con_x_st
        Graphe_distances[(-1, k)] = 0 if len(con_x_st) > len(pro_x_st) else (
                len(con_x_st) + (1 if len(con_x_st) < len(pro_x_st) else 0))

        st_y_cofactor = st[1] - alt_y
        pro_st_y = {i_ for i_ in range(len(st_y_cofactor)) if st_y_cofactor[i_] == 1}
        con_st_y = {j_ for j_ in range(len(st_y_cofactor)) if st_y_cofactor[j_] == -1}
        Pros_Dict[(k, -2)] = pro_st_y
        Cons_Dict[(k, -2)] = con_st_y
        Graphe_distances[(k, -2)] = 0 if len(con_st_y) > len(pro_st_y) else (
                len(con_st_y) + (1 if len(con_st_y) < len(pro_st_y) else 0))

    # liaison directe x->y
    x_y_cofactor = alt_x - alt_y
    pro_x_y = {i_ for i_ in range(len(x_y_cofactor)) if x_y_cofactor[i_] == 1}
    con_x_y = {j_ for j_ in range(len(x_y_cofactor)) if x_y_cofactor[j_] == -1}
    Pros_Dict[(-1, -2)] = pro_x_y
    Cons_Dict[(-1, -2)] = con_x_y
    Graphe_distances[(-1, -2)] = 0 if len(con_x_y) > len(pro_x_y) else (
            len(con_x_y) + (1 if len(con_x_y) < len(pro_x_y) else 0))

    model_gurobi = Model("Calcul PCC for explanations")
    model_gurobi.setParam('OutputFlag', False)
    nb_criteria_considered_alias_m = len(to_explain_statement[0])
    model_gurobi.Params.FeasibilityTol = CONSTRAINTSFEASIBILITYTOL
    Var_omega = [model_gurobi.addVar(vtype=GRB.CONTINUOUS, lb=CONSTRAINTSFEASIBILITYTOL * 10, name="w_{}".format(i_))
                 for
                 i_ in
                 range(nb_criteria_considered_alias_m)]

    Arcs_Liaisons_Var_Dict = {arc: model_gurobi.addVar(vtype=GRB.BINARY, name=f'arc_{arc}') for arc in Graphe_distances}

    # contrainte : exactement un arc sortant du sommet -1
    model_gurobi.addConstr(quicksum([var_b for arc, var_b in Arcs_Liaisons_Var_Dict.items() if arc[0] == -1]) == 1)

    # contraintes : continuité du chemin
    for k in range(len(PI_statements_List)):
        model_gurobi.addConstr(
            quicksum([var_b for arc, var_b in Arcs_Liaisons_Var_Dict.items() if arc[0] == k]) == quicksum(
                [var_b for arc, var_b in Arcs_Liaisons_Var_Dict.items() if arc[1] == k]))

    B_i_j_VarDict = {
        arc: {(i_, j_): model_gurobi.addVar(vtype=GRB.BINARY, name=f'b_{i_}_{j_}_{arc}') for i_ in Pros_Dict[arc] for j_
              in Cons_Dict[arc]} for arc in Graphe_distances}

    for arc in Graphe_distances:
        for i_ in Pros_Dict[arc]:  # pro used at most once
            model_gurobi.addConstr(
                quicksum([var_concerned for (i_b, j_b), var_concerned in B_i_j_VarDict[arc].items() if i_ == i_b]) <=
                Arcs_Liaisons_Var_Dict[arc])

        for j_ in Cons_Dict[arc]:  # cons used once
            model_gurobi.addConstr(
                quicksum([var_concerned for (i_b, j_b), var_concerned in B_i_j_VarDict[arc].items() if j_ == j_b]) ==
                Arcs_Liaisons_Var_Dict[arc])

        for (i_, j_), var_concerned in B_i_j_VarDict[arc].items():
            model_gurobi.addConstr(
                Var_omega[i_] - Var_omega[j_] >= (1 + CONSTRAINTSFEASIBILITYTOL * 10) * var_concerned - 1)

    for elmt in PI_cofactors_List:
        model_gurobi.addConstr(quicksum(np.array(Var_omega) * elmt) >= CONSTRAINTSFEASIBILITYTOL)

    model_gurobi.addConstr(quicksum(Var_omega) == 1.)

    model_gurobi.params.PoolSolutions = poolSolutions
    model_gurobi.params.PoolSearchMode = 2
    model_gurobi.params.PoolGap = len(PI_statements_List)

    model_gurobi.update()
    model_gurobi.setObjective(
        quicksum([(1 + Graphe_distances[arc]) * var_b for arc, var_b in Arcs_Liaisons_Var_Dict.items()]), GRB.MINIMIZE)
    model_gurobi.optimize()

    if model_gurobi.status == GRB.OPTIMAL:
        All_explanations = list()
        for sol_number in range(model_gurobi.SolCount):
            model_gurobi.params.SolutionNumber = sol_number
            List_Arcs_PCC = [-1]
            courant = List_Arcs_PCC[-1]
            while courant != -2:
                suiv_courant = -2
                for k in range(len(PI_statements_List)):
                    if (courant, k) in Arcs_Liaisons_Var_Dict and round(
                            Arcs_Liaisons_Var_Dict[(courant, k)].Xn) == 1:
                        suiv_courant = k
                        List_Arcs_PCC.append(suiv_courant)
                courant = suiv_courant
            List_Arcs_PCC.append(-2)

            explList = list()
            for j in range(1, len(List_Arcs_PCC)):
                explList.append([(i_, j_) for (i_, j_), var_concerned in
                                 B_i_j_VarDict[(List_Arcs_PCC[j - 1], List_Arcs_PCC[j])].items() if
                                 round(var_concerned.Xn) == 1 and Var_omega[i_].Xn - Var_omega[j_].Xn >= 0])
                if len(explList[-1]) < len(Cons_Dict[(List_Arcs_PCC[j - 1], List_Arcs_PCC[j])]):
                    break
                explList.append(List_Arcs_PCC[j])
            else:
                All_explanations.append(explList[:-1])

        return True, All_explanations

    return False, None





