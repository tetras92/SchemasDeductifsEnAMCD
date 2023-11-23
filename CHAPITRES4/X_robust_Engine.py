import numpy as np

from Reference_Functions_Robuste import *


def compute_robust_explanations(Statements_Dict, PI_statements_List, PI_cofactors_List, evolutive_context=True):
    PairwiseXplainedWithDetailsDict = dict()
    Xplained_statements = list()
    Xplained_cofactors = list()

    step = 0
    stabilised = False
    Schemes_used = set()
    BDB = {"[DoPa]": dict(), "[CePa]": dict(), "[Tr]": dict(), "[TrE]": dict(), "[Cov]": dict(), "[CovD]": dict(),
           "[1-Cg]": dict(), "[1-CgD]": dict(), "[2-Cg]": dict(), "[2-CgD]": dict()}
    while not stabilised:
        # print(len(PairwiseXplainedWithDetailsDict))
        stabilised = evolutive_context
        TmpXplained_statements = list()
        TmpXplained_cofactors = list()

        for (k1, k2), statementk1k2 in Statements_Dict.items():
            if (k1, k2) in PairwiseXplainedWithDetailsDict:
                continue

            cofactork1k2 = statementk1k2[0] - statementk1k2[1]

            is_pareto, d_pareto = doPa_robust(cofactork1k2, step=step, BDB=BDB)
            if is_pareto:
                stabilised = not evolutive_context
                PairwiseXplainedWithDetailsDict[(k1, k2)] = d_pareto
                TmpXplained_statements.append(statementk1k2)
                TmpXplained_cofactors.append(cofactork1k2)
                # Schemes_used.add("[DoPa]")
                continue

            is_ct, d_ct = cePa_robust(statementk1k2, PI_statements_List + Xplained_statements,
                                      PI_cofactors_List + Xplained_cofactors, step=step, BDB=BDB)
            if is_ct:
                stabilised = not evolutive_context
                PairwiseXplainedWithDetailsDict[(k1, k2)] = d_ct
                TmpXplained_statements.append(statementk1k2)
                TmpXplained_cofactors.append(cofactork1k2)
                # Schemes_used.add("[CePa]")
                continue

            is_tr, d_tr = tr_robust(statementk1k2, PI_statements_List + Xplained_statements, step=step, BDB=BDB)
            if is_tr:
                stabilised = not evolutive_context
                PairwiseXplainedWithDetailsDict[(k1, k2)] = d_tr
                TmpXplained_statements.append(statementk1k2)
                TmpXplained_cofactors.append(cofactork1k2)
                # Schemes_used.add("[Tr]")
                continue

            is_trp, d_trp = trE_robust(statementk1k2, PI_statements_List + Xplained_statements, step=step, BDB=BDB)
            if is_trp:
                stabilised = not evolutive_context
                PairwiseXplainedWithDetailsDict[(k1, k2)] = d_trp
                TmpXplained_statements.append(statementk1k2)
                TmpXplained_cofactors.append(cofactork1k2)
                # Schemes_used.add("[TrE]")
                continue

            is_kb, d_kb = cov0_robust(cofactork1k2, PI_cofactors_List + Xplained_cofactors, step=step, BDB=BDB)
            if is_kb:
                stabilised = not evolutive_context
                PairwiseXplainedWithDetailsDict[(k1, k2)] = d_kb
                TmpXplained_statements.append(statementk1k2)
                TmpXplained_cofactors.append(cofactork1k2)
                # Schemes_used.add("[Cov]")
                continue

            is_gns, d_gns = covD0_robust(statementk1k2, PI_statements_List + Xplained_statements,
                                         PI_cofactors_List + Xplained_cofactors, step=step, BDB=BDB)
            if is_gns:
                stabilised = not evolutive_context
                PairwiseXplainedWithDetailsDict[(k1, k2)] = d_gns
                TmpXplained_statements.append(statementk1k2)
                TmpXplained_cofactors.append(cofactork1k2)
                # Schemes_used.add("[CovD]")
                continue

            is_cg_all, d_cg_all = un_Cg_robust(statementk1k2, PI_statements_List + Xplained_statements,
                                               PI_cofactors_List + Xplained_cofactors, step=step, BDB=BDB)
            if is_cg_all:
                stabilised = not evolutive_context
                PairwiseXplainedWithDetailsDict[(k1, k2)] = d_cg_all
                TmpXplained_statements.append(statementk1k2)
                TmpXplained_cofactors.append(cofactork1k2)
                # Schemes_used.add("[1-Cg]")
                continue

            is_m_cg_all, d_m_cg_all = un_Cg_D_robust(statementk1k2, PI_statements_List + Xplained_statements,
                                                     PI_cofactors_List + Xplained_cofactors, step=step, BDB=BDB)
            if is_m_cg_all:
                stabilised = not evolutive_context
                PairwiseXplainedWithDetailsDict[(k1, k2)] = d_m_cg_all
                TmpXplained_statements.append(statementk1k2)
                TmpXplained_cofactors.append(cofactork1k2)
                # Schemes_used.add("[1-CgD]")
                continue

            is_reduced_m_cg_all, d_reduced_m_cg_all = deux_Cg_robust(statementk1k2,
                                                                     PI_statements_List + Xplained_statements,
                                                                     PI_cofactors_List + Xplained_cofactors, step=step, BDB=BDB)
            if is_reduced_m_cg_all:
                stabilised = not evolutive_context
                PairwiseXplainedWithDetailsDict[(k1, k2)] = d_reduced_m_cg_all
                TmpXplained_statements.append(statementk1k2)
                TmpXplained_cofactors.append(cofactork1k2)
                # Schemes_used.add("[2-Cg]")
                continue

            is_2_cg_d, d_2_cg_d = deux_Cg_D_robust(statementk1k2, PI_statements_List + Xplained_statements,
                                                   PI_cofactors_List + Xplained_cofactors, step=step, BDB=BDB)
            if is_2_cg_d:
                stabilised = not evolutive_context
                PairwiseXplainedWithDetailsDict[(k1, k2)] = d_2_cg_d
                TmpXplained_statements.append(statementk1k2)
                TmpXplained_cofactors.append(cofactork1k2)
                Schemes_used.add("[2-CgD]")
                continue

        Xplained_statements += TmpXplained_statements
        Xplained_cofactors += TmpXplained_cofactors
        if not stabilised:
            To_remove = [(cat, elem) for cat, DCat in BDB.items() for elem, valueslist in DCat.items() if
                         not valueslist[0]]
            for cat, elem in To_remove:
                del BDB[cat][elem]
        step += 1
    return PairwiseXplainedWithDetailsDict


