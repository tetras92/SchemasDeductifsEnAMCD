import numpy as np
from Utils import check_int_List_1, check_int_List_2, integer_to_bin


def generate_instance(m, size_A=10, size_AR=3):
    """m : nombre de critères; size_A: nombre d'alternatives du problème; size_AR : nombre d'alternatives de référence"""

    L = [round(float(np.random.rand()), 8) for _ in range(m)]                 # 18/09/22 : S'ASSURER QUE LE CONSTRAINTSFEASIBILITYTOL EST TOUJOURS À 10-9
    L = sorted(L) + [1.]
    omega_score = np.array([L[i] - L[i - 1] for i in range(1, len(L))])
    # print(omega_score)
    setA_as_integers = tuple(*np.random.choice(range(1, 2 ** m - 1), size=(1, size_A), replace=False))
    while not (check_int_List_2(setA_as_integers, m) and check_int_List_1(setA_as_integers, m)):
        setA_as_integers = tuple(*np.random.choice(range(1, 2 ** m - 1), size=(1, size_A), replace=False))

    setA_as_binary_np_arrays = [np.array([int(binary_digit) for binary_digit in integer_to_bin(int(alt_int_val), m)])
                                for alt_int_val in setA_as_integers]
    setA_as_binary_np_arrays = list(
        sorted(setA_as_binary_np_arrays, key=lambda arr: sum(arr * omega_score), reverse=True))

    setAR_as_integers = tuple(*np.random.choice(range(1, 2 ** m - 1), size=(1, size_AR), replace=False))
    while not (check_int_List_2(setAR_as_integers, m) and check_int_List_1(setAR_as_integers, m)):
        setAR_as_integers = tuple(*np.random.choice(range(1, 2 ** m - 1), size=(1, size_AR), replace=False))

    setAR_as_binary_np_arrays = [np.array([int(binary_digit) for binary_digit in integer_to_bin(int(alt_int_val), m)])
                                 for alt_int_val in setAR_as_integers]
    setAR_as_binary_np_arrays = list(
        sorted(setAR_as_binary_np_arrays, key=lambda arr: sum(arr * omega_score), reverse=True))

    # print("w : ", omega_score)
    # print("A : ")
    # print(*setA_as_binary_np_arrays, sep="\n")
    # print("AR : ")
    # print(*setAR_as_binary_np_arrays, sep="\n")

    return omega_score, setA_as_binary_np_arrays, setAR_as_binary_np_arrays

# generate_instance(6)

