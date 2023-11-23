from gurobipy import *


def decompose(proSet, conSet, W):
    model = Model("PLNE for Delta 11 decomposition")
    model.setParam('OutputFlag', False)
    #-- Variables
    Sij = {(i, j): model.addVar(vtype=GRB.BINARY, name=f's_{i}_{j}') for i in proSet for j in conSet}


    # Couverture de tout critere con
    for j in conSet:
        model.addConstr(quicksum([Sij[(i_, j)] for i_ in proSet]) == 1)

    # Chaque critere pro utilise au plus une fois
    for i in proSet:
        model.addConstr(quicksum([Sij[(i, j_)] for j_ in conSet]) <= 1)

    # Compatibilite swap et jeu de poids
    for i in proSet:
        for j in conSet:
            model.addConstr(W[j]*Sij[(i, j)] <= W[i])

    model.update()
    model.optimize()

    if model.status == GRB.OPTIMAL:
        chainons_arguments_list = list()
        dominance_swap = {i: ({i}, set()) for i in proSet}
        for j in conSet:
            for i in proSet:
                if Sij[(i, j)].x == 1:
                    chainons_arguments_list.append(({i}, {j}))
                    del dominance_swap[i]

        return True, chainons_arguments_list + list(dominance_swap.values())
    else:
        return False, list()


