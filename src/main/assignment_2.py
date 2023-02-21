#Author: Samuel Yunker
import numpy as np


def LagrangianCoefficient(n:int,k:int,x:float,points:np.ndarray):
    """
    Returns the L_{n,k}(x) for a given set of points
    parameters:
        - n(int): The total number of points to use
            - when n = len(points) + 1 all points will be used
        - k(int): the lagrangian coeficient index
        - x(float): the current point to get the coeficient of
        - points(np.ndarray): 1D array of sample points given for the Lagrangian
    """
    cumProduct = 1
    for i in range(n):
        if(i==k):
            continue
        else:
            currentProduct = (x-points[i]) / (points[k] - points[i])
            cumProduct *= currentProduct

    return cumProduct

def LagrangianPolynomial(n:int,x:float,points_x:np.ndarray,points_y:np.ndarray):
    """
    Returns the Lagrangian Polynomial P_{n}(x) for a given set of points.
    parameters:
        - n(int): The total number of points to use
            - when n = len(points) + 1 all points will be used
        - x(float): the current point
        - points_x(np.ndarray): 1D array of sample x points given for the Lagrangian
        - points_y(np.ndarray): 1D array of sample y points given for the Lagrangian
    """

    cumSum = 0
    for k in range(n):
        currentSum = points_y[k] * LagrangianCoefficient(n,k,x,points_x)
        cumSum += currentSum
    
    return cumSum


def nevillesMethod(points_x:np.ndarray, points_y:np.ndarray, x:float):

    n = min(len(points_x),len(points_y))
    nevillesTable = np.zeros((n, n))

    # fill in value (just the y values because we already have x set)
    for i in range(n):
        nevillesTable[i,0] = points_y[i]

    for j in range(1,n):
        for i in range(j,n):
            p0 = nevillesTable[i-1,j-1]
            x0 = points_x[i-j]

            p1 = nevillesTable[i,j-1]
            x1 = points_x[i]

            coef = 1 / (x1 - x0)
            mult1 = (x-x0)*p1
            mult2 = (x-x1)*p0
            tableValue = coef*(mult1 - mult2)
            nevillesTable[i,j] = tableValue
    #print(nevillesTable)
    return nevillesTable[-1,-1]


    
    
    return 0

def buildDDTable(points_x:np.ndarray, points_y:np.ndarray):
    """
    Creates and returns a divided difference table for a given set of points.
    parameters:
        - points_x(np.ndarray): 1D array of sample x points
        - points_y(np.ndarray): 1D array of sample y points
    """

    n = min(len(points_x),len(points_y))
    dividedDifferenceTable = np.zeros((n, n))

    # fill in value (just the y values because we already have x set)
    for i in range(n):
        dividedDifferenceTable[i,0] = points_y[i]

    for j in range(1,n):
        for i in range(j,n):
            f0 = dividedDifferenceTable[i-1,j-1]
            x0 = points_x[i-j]

            f1 = dividedDifferenceTable[i,j-1]
            x1 = points_x[i]

            
            numerator = f1 - f0
            denominator = x1-x0
            tableValue = numerator/denominator
            dividedDifferenceTable[i,j] = tableValue

    #print(dividedDifferenceTable)
    return dividedDifferenceTable

def newtonForward(dividedDifferenceTable:np.ndarray, points_x:np.ndarray, x:float,n:int):

    polynomial = dividedDifferenceTable[0,0]
    reoccuring_x_span = 1
    for i in range(1,n+1):
        coef = dividedDifferenceTable[i,i]

        x_root = x - points_x[i-1]
        reoccuring_x_span *= x_root

        polynomial += coef*reoccuring_x_span

    return polynomial

def buildHermiteDDTable(points_x:np.ndarray, points_y:np.ndarray, points_y_prime:np.ndarray):
    """
    Creates and returns a divided difference table for the hermite interpolation given set of points and the derivative at those points.
    parameters:
        - points_x(np.ndarray): 1D array of sample x points
        - points_y(np.ndarray): 1D array of sample y points
        - points_y_prime(np.ndarray): 1D array of derivatives of the function at each x,y point
    """

    n = 2*min(len(points_x),len(points_y))
    hermiteTable = np.zeros((n, n+1))

    # fill in value (just the y values because we already have x set)
    for i in range(n):
        hermiteTable[i,0] = points_x[i//2]
        hermiteTable[i,1] = points_y[i//2]

    #print(hermiteTable)
    for j in range(2,n+1):
        for i in range(j-1,n):
            if(j == 2 and i % 2 == 1):
                hermiteTable[i,2] = points_y_prime[i//2]
            else:
                f0 = hermiteTable[i-1,j-1]
                x0 = hermiteTable[i-(j-1),0]

                f1 = hermiteTable[i,j-1]
                x1 = hermiteTable[i,0]

                
                numerator = f1 - f0
                denominator = x1-x0
                tableValue = numerator/denominator
                hermiteTable[i,j] = tableValue

    #print(dividedDifferenceTable)
    return hermiteTable

def hermitePolynomialAproximation(points_x:np.ndarray, points_y:np.ndarray, points_y_prime:np.ndarray, x:float, n:int):
    DDtable = buildHermiteDDTable(points_x,points_y,points_y_prime)
    
    polynomial = DDtable[0,1]
    reoccuring_x_span = 1
    #print(f"p_{n}({x}) = {polynomial} ",end='')
    for i in range(1,n+1):
        coef = DDtable[i,i+1]

        x_root = x - DDtable[i-1,0]
        reoccuring_x_span *= x_root
        #print(f"+ {coef}*{reoccuring_x_span} ",end='')

        polynomial += coef*reoccuring_x_span
    #print(f"\n\t={polynomial}")
    

    return polynomial


def cubicSpline():
    return 0


def question1():

    points_x = np.array([3.6,3.8,3.9])
    points_y = np.array([1.675,1.436,1.318])
    x = 3.7

    result = nevillesMethod(points_x,points_y,x)
    print(result)

def question2():

    points_x = np.array([7.2,7.4,7.5,7.6])
    points_y = np.array([23.5492,25.3913,26.8224,27.4589])

    #points_x = np.array([1,1.3,1.6,1.9,2.2])
    #points_y = np.array([0.7651977,0.6200860,0.4554022,0.2818186,0.1103623])
    x = 7.3

    ddTable = buildDDTable(points_x,points_y)
    polynomial_degree1 = newtonForward(ddTable,points_x,x,1)
    polynomial_degree2 = newtonForward(ddTable,points_x,x,2)
    polynomial_degree3 = newtonForward(ddTable,points_x,x,3)
    
    print(ddTable)
    print()
    print(polynomial_degree1)
    print(polynomial_degree2)
    print(polynomial_degree3)

def question4():
    # points_x = np.array([1.3,1.6,1.9])
    # points_y = np.array([0.6200860,0.4554022,0.2818186])
    # points_y_prime = np.array([-0.5220232, -0.5698959, -0.5811571])
    # x = 1.5
    points_x = np.array([3.6,3.8,3.9])
    points_y = np.array([1.675,1.436,1.318])
    points_y_prime = np.array([-1.195,-1.188,-1.182])

    DDtable = buildHermiteDDTable(points_x,points_y,points_y_prime)
    #np.set_printoptions(precision=5)
    print(DDtable[:,:5])

    #result = hermitePolynomialAproximation(points_x,points_y,points_y_prime,x,5)
    #print(result)

if __name__ == "__main__":
    #question1()
    #question2()
    question4()
    