import csv
import itertools as iter_

from Reference_Functions_lv import *

m_value = 4
directory_CBTO = f'/home/manuel239/PycharmProjects/SchemasDeductifsEnAMCD/CBTO{m_value}/'


def est_compatible_avec(omega_array, PI_statements_arrays):
    for statement_array in PI_statements_arrays:
        X_array, Y_array = statement_array
        if sum(omega_array * X_array) < sum(omega_array * Y_array):
            return False
    return True


bc_ad = (np.array([0, 1, 1, 0]), np.array([1, 0, 0, 1]))
c_ab = (np.array([0, 0, 1, 0]), np.array([1, 1, 0, 0]))
b_a = (np.array([0, 1, 0, 0]), np.array([1, 0, 0, 0]))
d_ac = (np.array([0, 0, 0, 1]), np.array([1, 0, 1, 0]))
ac_b = (np.array([1, 0, 1, 0]), np.array([0, 1, 0, 0]))
# c_ab = (np.array([0, 0, 1, 0]), np.array([1, 1, 0, 0]))
# d_b  = (np.array([0, 0, 0, 1]), np.array([0, 1, 0, 0]))
b_ac = (np.array([0, 1, 0, 0]), np.array([1, 0, 1, 0]))
d_ab = (np.array([0, 0, 0, 1]), np.array([1, 1, 0, 0]))
PI = [bc_ad, d_ac]#, d_ab]


def rend_vainqueur_sachant(omega_array, Dominance_arrays):
    return est_compatible_avec(omega_array, Dominance_arrays)


bd_ac = (np.array([0, 1, 0, 1]), np.array([1, 0, 1, 0]))
d_ab =  (np.array([0, 0, 0, 1]), np.array([1, 1, 0, 0]))
cd_a =  (np.array([0, 0, 1, 1]), np.array([1, 0, 0, 0]))

O1 = [bd_ac, d_ab, cd_a]

ac_bd = (np.array([1, 0, 1, 0]), np.array([0, 1, 0, 1]))
ac_b  = (np.array([1, 0, 1, 0]), np.array([0, 1, 0, 0]))
acd_b = (np.array([1, 1, 0, 1]), np.array([0, 1, 0, 0]))

O2 = [ac_bd, ac_b, acd_b]

ab_d = (np.array([1, 1, 0, 0]), np.array([0, 0, 0, 1]))
b_ac = (np.array([0, 1, 0, 0]), np.array([1, 0, 1, 0]))

O3 = [ab_d, b_ac]

if __name__ == "__main__":
    CPT = 0
    TOTAL = 0
    for id_model in range(1, len(os.listdir(directory_CBTO)) + 1):
        omega_file = f'model{id_model}.csv'
        with open(directory_CBTO + omega_file) as utilityFile:
            reader = csv.DictReader(utilityFile)
            w_list = list()
            for row in reader:
                for criterion in reader.fieldnames:
                    w_list.append(float(row[criterion]))

        for w_list_perm in iter_.permutations(w_list):
            omega_score_function = np.array(w_list_perm)
            omega_score_function_dict = {i: omega_score_function[i] for i in range(len(omega_score_function))}
            TOTAL += 1
            if est_compatible_avec(omega_score_function, PI):
                print(id_model, {chr(ord('a') + i): val for i, val in omega_score_function_dict.items()})
                CPT += 1
                if rend_vainqueur_sachant(omega_score_function, O1):
                    print("=>O1")
                if rend_vainqueur_sachant(omega_score_function, O2):
                    print("=>O2")
                if rend_vainqueur_sachant(omega_score_function, O3):
                    print("=>O3")
    print('NB total', CPT, "/", TOTAL)

# PI = [bc_ad, d_ac]
# 1 {'a': 1.0, 'b': 8.0, 'c': 2.0, 'd': 4.0}
# =>O3
# 1 {'a': 2.0, 'b': 8.0, 'c': 1.0, 'd': 4.0}
# =>O3
# 2 {'a': 2.0, 'b': 10.0, 'c': 3.0, 'd': 6.0}
# =>O3
# 2 {'a': 3.0, 'b': 10.0, 'c': 2.0, 'd': 6.0}
# =>O3
# 3 {'a': 2.0, 'b': 10.0, 'c': 4.0, 'd': 7.0}
# =>O3
# 3 {'a': 4.0, 'b': 10.0, 'c': 2.0, 'd': 7.0}
# =>O3
# 4 {'a': 1.0, 'b': 4.0, 'c': 6.0, 'd': 8.0}
# =>O1
# 4 {'a': 1.0, 'b': 6.0, 'c': 4.0, 'd': 8.0}
# =>O1
# 4 {'a': 1.0, 'b': 8.0, 'c': 4.0, 'd': 6.0}
# =>O3
# 5 {'a': 3.0, 'b': 10.0, 'c': 4.0, 'd': 8.0}
# =>O3
# 5 {'a': 4.0, 'b': 10.0, 'c': 3.0, 'd': 8.0}
# =>O3
# 6 {'a': 2.0, 'b': 7.0, 'c': 4.0, 'd': 8.0}
# =>O3
# 6 {'a': 2.0, 'b': 8.0, 'c': 4.0, 'd': 7.0}
# =>O3
# 10 {'a': 2.0, 'b': 6.0, 'c': 7.0, 'd': 10.0}
# =>O1
# 10 {'a': 2.0, 'b': 7.0, 'c': 6.0, 'd': 10.0}
# =>O1
# 12 {'a': 3.0, 'b': 8.0, 'c': 6.0, 'd': 10.0}
# NB total 16 / 336
#
# Process finished with exit code 0
