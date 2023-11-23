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


