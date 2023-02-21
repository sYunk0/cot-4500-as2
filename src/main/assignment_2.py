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


def nevillesMethod():
    return 0

def newtonForward():
    return 0

def dividedDiference():
    return 0

def cubicSpline():
    return 0


def question1():
    result = nevillesMethod()
    print(result)


if __name__ == "__main__":
    question1()