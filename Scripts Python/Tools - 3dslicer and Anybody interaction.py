# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 12:41:17 2024

Scripts for interactions between Anybody and 3d slicer

@author: Dan
"""

# %% Tools for converting codes between Anybody and python

def ExportMat2Anybody(transform_node_name: str):
    """
    Prints a 4x4 slicer rotation matrix from a Transform node to a code in the Anybody language for the rotation matrix and translation vector

    NodeName : str - name of the Transform node matrix

    Print
    ---------------
    Prints an AnyMat33 code for Anybody containing the rotation matrix
    Prints a AnyVec3 code for Anybody containing the translation vector

    """

    MatrixNode = slicer.util.getNode(transform_node_name)
    Matrix = vtk.vtkMatrix4x4()
    MatrixNode.GetMatrixTransformToWorld(Matrix)

    Matrix = slicer.util.arrayFromVTKMatrix(Matrix)

    Matrix = np.around(Matrix, decimals=6)

    Rotation = f'AnyMat33 Rotation = {{\n{{{Matrix[0,0]} , {Matrix[0,1]} , {Matrix[0,2]}}},\n{{{Matrix[1,0]} , {Matrix[1,1]} , {Matrix[1,2]}}},\n{{{Matrix[2,0]} , {Matrix[2,1]} , {Matrix[2,2]}}}\n\n}};'
    Translation = f'AnyVec3 Position = 0.001*{{{Matrix[0,3]} , {Matrix[1,3]} , {Matrix[2,3]}}};'

    print(Translation)
    print(Rotation)


def AnyMatrix2Array(string):
    """
    Converts a string containing the anybody definition of a matrix and converts it to an numpy.array
    """
    string = string.replace(" ", "")
    matrix = []
    if string[0:2] == "{{":
        Type = "Matrix"
        rows = string.split("},{")
        rows[0] = rows[0].replace("{{", "")
        rows[-1] = rows[-1].replace("}}", "")
    else:

        rows = string.split("},{")
        rows[0] = rows[0].replace("{", "")
        rows[-1] = rows[-1].replace("}", "")

        if len(rows[0].split(",")) == 1:
            Type = "Value"
        else:
            Type = "Vector"

    if Type == "Value":
        matrix = eval(rows[0])
    elif Type == "Vector":
        rows = rows[0].split(",")
        rows = [eval(i) for i in rows]
        matrix.append(rows)
        matrix = np.array(matrix)
    else:
        for row in rows:
            elements = row.split(",")
            elements = [eval(i) for i in elements]
            matrix.append(elements)
        matrix = np.array(matrix)

    return matrix


# %% Tools for 3d slicer

def TransformMatrix_4x4(RotMat=None, Vect=None):
    """
    Input : rotation matrix and translation vector to a 4x4 matrix for a 3d slicer transform
    """

    if RotMat is None:
        RotMat = np.eye(3)
    if Vect is None:
        Vect = np.zeros(3)

    Mat = np.eye(4)

    Mat[0:3, 0:3] = RotMat
    Mat[0:3, 3] = Vect
    return Mat


def Set_TransformMatrix_4x4(TransformName, RotMat=None, Vect=None):
    """
    Function to set the 4x4 transform matrix to a transform in 3d slicer from a rotation matrix and translation vector numpy arrays
    """

    TransformMat = TransformMatrix_4x4(RotMat, Vect)

    Transform = getNode(TransformName)
    Transform.SetMatrixTransformToParent(slicer.util.vtkMatrixFromArray(TransformMat))