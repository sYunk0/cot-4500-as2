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

def newtonForward():
    return 0

def dividedDiference():
    return 0

def cubicSpline():
    return 0


def question1():

    points_x = np.array([3.6,3.8,3.9])
    points_y = np.array([1.675,1.436,1.318])
    x = 3.7

    result = nevillesMethod(points_x,points_y,x)
    print(result)


if __name__ == "__main__":
    question1()
    