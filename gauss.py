eps = 1e-10    # погрешность

# Функция вывода матрицы на экран
def displayArray(info ,a):
    print(info)
    n = len(a)
    for i in range(n):
        line = ""
        for j in range(n + 1):
            line += str("%10.5f" % a[i][j]) + "      "
            if j == n - 1:
                line += "| "
        print(line)
    print("")


# Функция вывода матрицы на экран
def displayArray2(info ,a):
    print(info)
    n = len(a)
    for i in range(n):
        line = ""
        for j in range(n):
            line += str("%10.5f" % a[i][j]) + "      "
        print(line)
    print("")


# Функция вывод решения на экран
def displaySolution(x):
    print("Решение системы:")
    for i, val in enumerate(x):
        print("x%x" % (i + 1), " = %5.5f" % val)


# Функция вывода умножения матриц на экран
def displayVec(info, r):
    print(info)
    n = len(r)
    for i in range(n):
        print("%5.5f" % r[i], end="    ")
    print("")


# Нахождение максимально элемента в строке
def maxelement(a, col, count_swap):
    n = len(a)
    maxel = a[col][col]
    maxrow = col
    for i in range(col + 1, n):
        if maxel < abs(a[i][col]):
            maxel = abs(a[i][col])
            maxrow = i
    if col != maxrow:
        swap(a, maxrow, col)
        count_swap += 1
    return count_swap


# Проверка на ноль
def checkByZero(q):
    if abs(q) < eps:
       return 1


# Вычисление определителя
def determinant(a, count_swap):
    det = 1
    n = len(a)
    if count_swap % 2:
        count_swap = -1
    else:
        count_swap = 1
    for i in range(n):
        det *= a[i][i]
    det *= count_swap
    return det


# Перестановка строк
def swap(a, row_one, row_two=0):
    n = len(a)
    for i in range(n + 1):
        tmp = a[row_one][i]
        a[row_one][i] = a[row_two][i]
        a[row_two][i] = tmp



# Вычетание строк
def sub(a, row_one, row_two, mn=1):
    n = len(a)
    for i in range(n + 1):
        a[row_one][i] -= a[row_two][i] * mn
    return a


# Приведение к треугольной матрице
def triangle(a):
    n = len(a)
    count_swap = 0
    for j in range(n):
        count_swap = maxelement(a, j, count_swap)
        for i in range(j + 1, n):
            c = a[i][j] / a[j][j]
            sub(a, i, j, c)
    return count_swap


# Нахождение решение, методом обратного хода
def searchSolution(a):
    n = len(a)
    solution = [0 for i in range(n)]
    for i in range(n - 1, -1, -1):
        solution[i] = a[i][n] / a[i][i]
        for j in range(i - 1, -1, -1):
            a[j][n] -= a[j][i] * solution[i]
    return solution


# Метод Гаусса
def gauss(a):
    count_swap = triangle(a)
    det = determinant(a, count_swap)
    flag = checkByZero(det)
    if flag:
        print("\nМатрица вырожденная. Определитель равен нулю\n")
        exit(1)
    x = searchSolution(a)
    return x


# Умножение матриц
def matrix_mul(a, x):
    n = len(a)
    m = len(x)
    result = []
    for i in range(n):
        s = 0
        for j in range(m):
            s += x[j] * a[i][j]
        result.append(s)
    return result


# Умножение матриц 2
def matrix_mul_2(a, b):
    n = len(a)
    res = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                res[i][j] += a[i][k] * b[k][j]
    return res

# Выделение матрицы B
def getVectorB(a):
    n = len(a)
    vectorB = []
    for i in range(n):
        vectorB.append(a[i][n])
    return vectorB

def getVectorA(a):
    vectorA = a[:]
    return vectorA

# Нахождение невязки
def discrepancy(res, b):
    for i in range(len(b)):
        res[i] = b[i] - res[i]
    return res


# Обратная матрица
def inverse(a):
    n = len(a)
    inv = [[0 for i in range(len(a) + 1)] for j in range(len(a))]
    for i in range(n):
        for j in range(n):
            if i == j:
                a[j][n] = 1.0
            else:
                a[j][n] = 0.0
        q = copy.deepcopy(a)
        vec_sol = gauss(q)
        for k in range(len(vec_sol)):
            inv[k][i] = vec_sol[k]
    return inv


# Вычисление нормы
def norm(a):
    n = len(a)
    max = 0
    for i in range(n):
        sum = 0
        for j in range(n):
            sum += abs(a[i][j])
            if max < sum:
                max = sum
    return max


# Число обусловленности матрицы
def condition_number(norm_one, norm_two):
    cond = norm_one * norm_two
    return cond

if __name__ == "__main__":
    import numpy, copy

    arr = numpy.loadtxt('input.txt', float)
    clean_arr = copy.deepcopy(arr)


    displayArray("Начальная матрица", arr)


    # Метод Гаусса
    count_swap = triangle(arr)
    displayArray("Треугольная матрица", arr)
    det = determinant(arr, count_swap)
    flag = checkByZero(det)
    if flag:
        print("\nМатрица вырожденная. Определитель равен нулю\n")
        exit(1)
    print("Определитель: %5.5f" % det)
    x = searchSolution(arr)
    displaySolution(x)

    # Нахождение невязки
    b = getVectorB(clean_arr)
    displayVec("Вектор B:", b)
    res = matrix_mul(clean_arr, x)
    displayVec("Результат перемножения матрицы A на вектор решения X:", res)
    dis = discrepancy(res, b)
    displayVec("Невязка:", dis)

    # Нахождение обратной матрицы
    inv = inverse(clean_arr)
    displayArray2("Обратная матрица", inv)

    # Нахождение числа обусловленности матрицы
    a = getVectorA(clean_arr)
    displayArray2("Матрица A", a)
    norm_a = norm(clean_arr)
    print("Норма А: ", norm_a, end="")
    norm_inv = norm(inv)
    print("\nНорма обратной матрицы: ", norm_inv, end="")
    cond = condition_number(norm_a, norm_inv)
    print("\nЧисло обусловленности матрицы: ", cond, end="")


