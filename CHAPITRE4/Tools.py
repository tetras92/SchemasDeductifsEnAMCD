import numpy as np
from termcolor import colored

COLOR = "blue"
attribute_creator = lambda criterion, attribute_value : '{}:{}'.format(criterion,attribute_value)
colored_character = lambda c, o, color : colored(c, color) if c != o else c
colored_expression = lambda alternative1, alternative2 : ("".join([colored_character(alternative1[i], alternative2[i], COLOR) for i in range(len(alternative1))]),
                                                          "".join([colored_character(alternative2[i], alternative1[i], COLOR) for i in range(len(alternative1))]))


tradeoff = lambda s, VarList: "[{bj_sup} -> {bj_inf}, {bi_inf} -> {bi_sup}]".format(bi_inf=VarList[s[0]][s[2]][0],
                                                                                     bi_sup=VarList[s[0]][s[3]][0],
                                                                                     bj_inf=VarList[s[1]][s[4]][0],
                                                                                     bj_sup=VarList[s[1]][s[5]][0])
symbol = lambda x : "+" if x == 1 else "-"
EPSILON = 0.000001                      #
CONSTRAINTSFEASIBILITYTOL = 1e-6  #0.000000001 # borne min dans Gurobi 1e-9
INTEGERFEASIBILITYTOL = 1e-6
ROUNDED_NUMBER_OF_DECIMALS = 4
EMPTYSET = "emptyset"
def covectorOfPairWiseInformationWith2Levels(coupleAlt):
    alt1, alt2 = coupleAlt
    cov = list()
    for level1i, level2i in list(zip(alt1.attributeLevelsList, alt2.attributeLevelsList)):
        if level1i < level2i:
            cov.append(-1)
        elif level1i > level2i:
            cov.append(1)
        else:
            cov.append(0)
    return np.array(cov)

def difficultyLevel(coupleAlt):
    alt1, alt2 = coupleAlt
    d_level = 0
    for level1i, level2i in list(zip(alt1.attributeLevelsList, alt2.attributeLevelsList)):
        if level1i != level2i:
            d_level += 1
    return d_level

