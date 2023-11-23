from gurobipy import *


def decompose(proSet, conSet, W):
    model = Model("PLNE for Delta 1m and Delta m1 decomposition")
    model.setParam('OutputFlag', False)
    #-- Variables
    B1m = {(i, j): model.addVar(vtype=GRB.BINARY, name=f'b_{i}_{j}_1m') for i in proSet for j in conSet}
    Bm1 = {(i, j): model.addVar(vtype=GRB.BINARY, name=f'b_{i}_{j}_m1') for i in proSet for j in conSet}

    # 1. tentative de modification (ce 16/02/22)
    for i in proSet:
        for j in conSet:
            model.addConstr(Bm1[(i, j)] + quicksum([B1m[(i_, j)] for i_ in proSet]) <= 1, name="C1 ({})".format((i, j)))

    # 2. Unicite contribution de i dans le Delta m1 world                 [Ce 28/11/21 jugé non indispensable à cause de 3. Unicite inter-worlds]
    # for i in proSet:
    #     model.addConstr(quicksum([Bm1[(i, j_)] for j_ in conSet]) <= 1, name="C2 ({})".format(i))

    # 3. Unicite inter-worlds
    for i in proSet:
        for j in conSet:
            model.addConstr(B1m[(i, j)] + quicksum([Bm1[(i, j_)] for j_ in conSet]) <= 1, name="C3 ({})".format((i, j)))

    # 4. Couverture Delta 1m
    for i in proSet:
        model.addConstr(quicksum([B1m[(i, j_)]*W[j_] for j_ in conSet]) <= W[i], name="C4 ({})".format(i))


    # 5. Couverture Delta m1
    for j in conSet:
        model.addConstr(quicksum([Bm1[(i_, j)]*W[i_] + B1m[(i_, j)]*W[i_] for i_ in proSet]) >= W[j], "C5 ({})".format(j))


    model.update()
    # model.setObjective(quicksum(B1m.values()), GRB.MAXIMIZE)
    model.optimize()

    if model.status == GRB.OPTIMAL:
        chainons_arguments_list = list()
        I1m = set()
        J1m = set()

        # print("\nDelta 1m")
        for i in proSet:
            if sum([int(B1m[(i, j_)].x) for j_ in conSet]) > 0:
                localJ1m = {j_ for j_ in conSet if int(B1m[(i, j_)].x) == 1}
                chainons_arguments_list.append(({i}, localJ1m))
                # print("{} -> {}".format(i, localJ1m))
                I1m.add(i)
                J1m = J1m | localJ1m

        # print("\nDelta m1")
        for j in conSet - J1m:
            if sum([int(Bm1[(i_, j)].x) for i_ in proSet]) > 0:
                localIm1 = {i_ for i_ in proSet if int(Bm1[(i_, j)].x) == 1}
                I1m = I1m | localIm1
                chainons_arguments_list.append((localIm1, {j}))

        if len(proSet - I1m) != 0:
            chainons_arguments_list.append((proSet - I1m, {}))

        return True, chainons_arguments_list
    else:
        return False, list()


