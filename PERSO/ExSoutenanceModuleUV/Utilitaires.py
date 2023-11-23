from gurobipy import *


def lafonction(nbUV, ContraintesExternes):
    mon_model = Model("Mon modele")
    E_UV_dict = {(i, j): mon_model.addVar(vtype=GRB.BINARY, name=f'e_uv{i}_{j}') for i in range(1, 5) for j in
                 range(1, nbUV + 1)}
    UV_M_dict = {(j, k): mon_model.addVar(vtype=GRB.BINARY, name=f'uv_m{j}_{k}') for j in range(1, nbUV + 1) for k in
                 range(1, 5)}

    E_E_dict = {(i1, i2, j): mon_model.addVar(vtype=GRB.BINARY, name=f'e_e{i1}_{i2}') for (i1, i2, _, _) in
                ContraintesExternes for j in range(1, nbUV + 1)}

    for (i1, i2, j), var_ee in E_E_dict.items():
        mon_model.addConstr(var_ee - E_UV_dict[(i1, j)] * E_UV_dict[(i2, j)] == 0)

    # 3 UV par modules
    for k in range(1, 5):
        mon_model.addConstr(quicksum([var_uv_m for jk, var_uv_m in UV_M_dict.items() if jk[1] == k]) >= 3)

    for j in range(1, nbUV + 1):
        mon_model.addConstr(quicksum([var_uv_m for jk, var_uv_m in UV_M_dict.items() if jk[0] == j]) >= 1)

    # contraintes externes
    for (i1, i2, k, aumoins2) in ContraintesExternes:
        premier_membre = LinExpr()
        for j in range(1, nbUV + 1):
            premier_membre += E_E_dict[(i1, i2, j)] * UV_M_dict[(j, k)]

        if aumoins2:
            mon_model.addConstr(premier_membre >= 2)
        else:
            mon_model.addConstr(premier_membre <= 1)
    mon_model.setParam('OutputFlag', False)
    mon_model.update()
    mon_model.optimize()

    return mon_model.status == GRB.OPTIMAL, {f'E{i}_UV{j}' for (i, j), euv in E_UV_dict.items() if euv.x == 1}, \
           {f'UV{j}_M{k}' for (j, k), jkv in UV_M_dict.items() if jkv.x == 1}




nbUV = 5
print(*lafonction(nbUV, [(2, 3, 1, True), (2, 4, 1, True), (1, 3, 1, True), (1, 4, 1, True), (3, 4, 1, False),
                        (1, 3, 2, True), (3, 4, 2, True), (1, 2, 2, True), (2, 4, 2, True), (1, 4, 2, False),
                        (1, 2, 3, True), (2, 3, 3, True), (1, 4, 3, True), (3, 4, 3, True), (1, 3, 3, False),
                        (1, 2, 4, True), (1, 3, 4, True), (2, 4, 4, True), (3, 4, 4, True), (2, 3, 4, False)]), sep="\n")


#=====================#
mon_model = Model("Mon modele")
Notes = {(i, j): mon_model.addVar(vtype=GRB.INTEGER, name=f'e_uv{i}_{j}', lb=11, ub=19) for i in range(1, 5) for j in range(1, 6)}
Coef = {j: mon_model.addVar(vtype=GRB.CONTINUOUS, name=f'coef{j}', lb=0.1, ub=0.9) for j in range(1, 6)}

# suivant M1 : [1, 2, 3]
M_UV = {1: {1, 2, 3}, 2: {1, 3, 4}, 3: {2, 3, 4}, 4: {2, 3, 5}}
E_UV = {1: {1, 2, 3, 5}, 2: {1, 2, 3, 4}, 3: {1, 3, 4, 5}, 4: {2, 3, 4, 5}}

Do = {1: [(2, 4), (2, 3), (3, 1), (4, 1)],
      2: [(3, 4), (3, 1), (4, 2), (1, 2)],
      3: [(2, 1), (2, 3), (1, 4), (3, 4)],
      4: [(1, 2), (1, 3), (2, 4), (3, 4)]
      }

for m, UVL in M_UV.items():
    mon_model.addConstr(quicksum([Coef[a] for a in UVL]) == 1.0)

for k, L in Do.items():
    UV_k = M_UV[k]
    for x, y in L:
        en_commun = E_UV[x] & E_UV[y]
        assert len(en_commun) >= 2
        mon_model.addConstr(quicksum(Notes[(x, a)] for a in en_commun) - quicksum(Notes[(y, a)] for a in en_commun) >= 0)
#
mon_model.setParam('OutputFlag', False)
mon_model.update()
mon_model.optimize()

if mon_model.status == GRB.OPTIMAL:
    print("Found !!!")
    for (i, j), nx in Notes.items():
        print((i, j), nx.x)
