from gurobipy import *
local_epsilon = 0.0001


def recommendation_after_deformation_L3(AlternativesSubsetsList, W_dict):
    SortedAlternativesSubsetsList = sorted(AlternativesSubsetsList,
                                           key=lambda alt: sum([W_dict[criterion] for criterion in alt]), reverse=True)

    model = Model("Deformation of DM weight -- Delta 1m-m1 -- L3 - with constraint above criteria order")
    DeviationsSigma = {k: (model.addVar(vtype=GRB.CONTINUOUS, lb=0, name=f'sig+_{k}'),
                           model.addVar(vtype=GRB.CONTINUOUS, lb=0, name=f'sig-_{k}'))
                       for k in W_dict}  # [0] = +; [1] = -

    model.setParam('OutputFlag', False)
    sum_w = 1. * sum(W_dict.values())
    W_dict = {k: vk / sum_w for k, vk in W_dict.items()}

    # ----------------------- CONSTRAINTS : SAME ORDER OVER CRITERIA ------------------------- #
    DecreasingOrderedCriteria = list(sorted(W_dict.keys(), key=lambda x: W_dict[x], reverse=True))
    for idk in range(1, len(DecreasingOrderedCriteria)):
        model.addConstr(W_dict[DecreasingOrderedCriteria[idk]] + DeviationsSigma[DecreasingOrderedCriteria[idk]][0] - DeviationsSigma[DecreasingOrderedCriteria[idk]][1] <=
                        W_dict[DecreasingOrderedCriteria[idk-1]] + DeviationsSigma[DecreasingOrderedCriteria[idk-1]][0] - DeviationsSigma[DecreasingOrderedCriteria[idk-1]][1] - local_epsilon)
    # ---------------------------------------------------------------------------------------- #

    E_ijVarDictPlus_I_1m = dict()
    E_ijVarDictMoins_I_1m = dict()
    E_ijVarDictPlus_J_1m = dict()
    E_ijVarDictMoins_J_1m = dict()

    E_ijVarDictPlus_m1 = dict()
    E_ijVarDictMoins_m1 = dict()

    other_alternative_related_swapVar_1m = dict()
    other_alternative_related_swapVar_m1 = dict()

    for v in range(1, len(SortedAlternativesSubsetsList)):
        oth_alt = SortedAlternativesSubsetsList[v]
        pros, cons = SortedAlternativesSubsetsList[0], oth_alt

        other_alternative_related_swapVar_1m[v] = {
            (i, j): model.addVar(vtype=GRB.BINARY, name=f's_{oth_alt}-{(i, j)}_1m') for i in pros for j in cons}
        other_alternative_related_swapVar_m1[v] = {
            (i, j): model.addVar(vtype=GRB.BINARY, name=f's_{oth_alt}-{(i, j)}_m1') for i in pros for j in cons}

        E_ijVarDictPlus_I_1m[v] = {(i, j): model.addVar(vtype=GRB.CONTINUOUS, lb=0, name=f'e+_{oth_alt}-{(i, j)}_I_1m')
                                   for i in pros for j in cons}
        E_ijVarDictMoins_I_1m[v] = {(i, j): model.addVar(vtype=GRB.CONTINUOUS, lb=0, name=f'e-_{oth_alt}-{(i, j)}_I_1m')
                                    for i in pros for j in cons}
        E_ijVarDictPlus_J_1m[v] = {(i, j): model.addVar(vtype=GRB.CONTINUOUS, lb=0, name=f'e+_{oth_alt}-{(i, j)}_J_1m')
                                   for i in pros for j in cons}
        E_ijVarDictMoins_J_1m[v] = {(i, j): model.addVar(vtype=GRB.CONTINUOUS, lb=0, name=f'e-_{oth_alt}-{(i, j)}_J_1m')
                                    for i in pros for j in cons}
        E_ijVarDictPlus_m1[v] = {(i, j): model.addVar(vtype=GRB.CONTINUOUS, lb=0, name=f'e+_{oth_alt}-{(i, j)}_m1') for
                                 i in pros for j in cons}
        E_ijVarDictMoins_m1[v] = {(i, j): model.addVar(vtype=GRB.CONTINUOUS, lb=0, name=f'e-_{oth_alt}-{(i, j)}_m1') for
                                  i in pros for j in cons}

        # -- Constraints
        for i in pros:
            model.addConstr(quicksum([other_alternative_related_swapVar_m1[v][(i, j_)] for j_ in cons]) <= 1)

        # - (b)
        for i in pros:
            for j in cons:
                model.addConstr(other_alternative_related_swapVar_1m[v][(i, j)] +
                                quicksum([other_alternative_related_swapVar_m1[v][(i, j_)] for j_ in cons]) <= 1)

        for i in pros:
            for j in cons:
                model.addConstr(other_alternative_related_swapVar_m1[v][(i, j)] + quicksum(
                    [other_alternative_related_swapVar_1m[v][(i_, j)] for i_ in pros]) <= 1)

        # - (d)
        for i in pros:
            w_i = W_dict[i]
            sig_i_plus = DeviationsSigma[i][0]
            sig_i_moins = DeviationsSigma[i][1]

            model.addConstr(w_i + sig_i_plus - sig_i_moins >= quicksum([other_alternative_related_swapVar_1m[v][
                                                                            (i, j_)] * W_dict[j_] +
                                                                        E_ijVarDictPlus_J_1m[v][(i, j_)] -
                                                                        E_ijVarDictMoins_J_1m[v][(i, j_)] for j_ in
                                                                        cons]))

        # - (e)
        for j in cons:
            w_j = W_dict[j]
            sig_j_plus = DeviationsSigma[j][0]
            sig_j_moins = DeviationsSigma[j][1]

            model.addConstr(quicksum([other_alternative_related_swapVar_1m[v][(i_, j)] * W_dict[i_] +
                                      other_alternative_related_swapVar_m1[v][(i_, j)] * W_dict[i_] +
                                      E_ijVarDictPlus_I_1m[v][(i_, j)] - E_ijVarDictMoins_I_1m[v][(i_, j)] +
                                      E_ijVarDictPlus_m1[v][(i_, j)] - E_ijVarDictMoins_m1[v][(i_, j)] for i_ in
                                      pros]) >= w_j + sig_j_plus - sig_j_moins)

        # linearisation
        for i in pros:
            sig_i_plus = DeviationsSigma[i][0]
            sig_i_moins = DeviationsSigma[i][1]
            for j in cons:
                sig_j_plus = DeviationsSigma[j][0]
                sig_j_moins = DeviationsSigma[j][1]

                model.addConstr(E_ijVarDictPlus_I_1m[v][(i, j)] <= other_alternative_related_swapVar_1m[v][(i, j)])
                model.addConstr(E_ijVarDictMoins_I_1m[v][(i, j)] <= other_alternative_related_swapVar_1m[v][(i, j)])
                model.addConstr(E_ijVarDictPlus_m1[v][(i, j)] <= other_alternative_related_swapVar_m1[v][(i, j)])
                model.addConstr(E_ijVarDictMoins_m1[v][(i, j)] <= other_alternative_related_swapVar_m1[v][(i, j)])

                model.addConstr(E_ijVarDictPlus_I_1m[v][(i, j)] <= sig_i_plus)
                model.addConstr(E_ijVarDictMoins_I_1m[v][(i, j)] <= sig_i_moins)
                model.addConstr(E_ijVarDictPlus_m1[v][(i, j)] <= sig_i_plus)
                model.addConstr(E_ijVarDictMoins_m1[v][(i, j)] <= sig_i_moins)

                model.addConstr(
                    E_ijVarDictPlus_I_1m[v][(i, j)] >= sig_i_plus - 1 + other_alternative_related_swapVar_1m[v][(i, j)])
                model.addConstr(
                    E_ijVarDictMoins_I_1m[v][(i, j)] >= sig_i_moins - 1 + other_alternative_related_swapVar_1m[v][
                        (i, j)])
                model.addConstr(
                    E_ijVarDictPlus_m1[v][(i, j)] >= sig_i_plus - 1 + other_alternative_related_swapVar_m1[v][(i, j)])
                model.addConstr(
                    E_ijVarDictMoins_m1[v][(i, j)] >= sig_i_moins - 1 + other_alternative_related_swapVar_m1[v][(i, j)])

                model.addConstr(E_ijVarDictPlus_J_1m[v][(i, j)] <= other_alternative_related_swapVar_1m[v][(i, j)])
                model.addConstr(E_ijVarDictMoins_J_1m[v][(i, j)] <= other_alternative_related_swapVar_1m[v][(i, j)])
                model.addConstr(E_ijVarDictPlus_J_1m[v][(i, j)] <= sig_j_plus)
                model.addConstr(E_ijVarDictMoins_J_1m[v][(i, j)] <= sig_j_moins)
                model.addConstr(
                    E_ijVarDictPlus_J_1m[v][(i, j)] >= sig_j_plus - 1 + other_alternative_related_swapVar_1m[v][(i, j)])
                model.addConstr(
                    E_ijVarDictMoins_J_1m[v][(i, j)] >= sig_j_moins - 1 + other_alternative_related_swapVar_1m[v][
                        (i, j)])

    model.addConstr(quicksum([sig_pair[0] for sig_pair in DeviationsSigma.values()]) == quicksum(
        [sig_pair[1] for sig_pair in DeviationsSigma.values()]))

    for k in W_dict:
        model.addConstr(DeviationsSigma[k][1] <= W_dict[k])

    model.update()
    model.setObjective(quicksum([sig_pair[0] for sig_pair in DeviationsSigma.values()]), GRB.MINIMIZE)
    model.optimize()
    if model.status == GRB.OPTIMAL:
        return True, float(model.objVal) == 0.0

    return False, None

