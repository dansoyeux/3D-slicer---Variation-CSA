# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 10:44:05 2023

@author: Dan

Scripts for CSA variation in 3d slicer
"""

import numpy as np
import math

"""
Glenoid inclination angle measure
"""

def get_angle_between_axes(axis_1, axis_2, name_angle=""):
    AxisNodeNames = [axis_1, axis_2]
    lineDirectionVectors = []
    for lineNodeName in AxisNodeNames:
        lineNode = slicer.util.getFirstNodeByClassByName("vtkMRMLMarkupsLineNode", lineNodeName)
        lineStartPos = np.zeros(3)
        lineEndPos = np.zeros(3)
        lineNode.GetNthControlPointPositionWorld(0, lineStartPos)
        lineNode.GetNthControlPointPositionWorld(1, lineEndPos)
        lineDirectionVector = (lineEndPos - lineStartPos) / np.linalg.norm(lineEndPos - lineStartPos)
        lineDirectionVectors.append(lineDirectionVector)
    angleRad = vtk.vtkMath.AngleBetweenVectors(lineDirectionVectors[0], lineDirectionVectors[1])
    angleDeg = vtk.vtkMath.DegreesFromRadians(angleRad)

    return angleDeg

# Print angles between slice nodes
def ShowAngle(axis_1, axis_2, name_angle=""):
    """
    prints the angle between 2 axes
    and prints in front the name of the angle
    """

    Angle = get_angle_between_axes(axis_1, axis_2, name_angle)

    print(name_angle + " = {0:0.3f}".format(Angle - 90) + "°")

    AxisNodeNames = [axis_1, axis_2]

    # Observe line node changes
    for lineNodeName in AxisNodeNames:
        lineNode = slicer.util.getFirstNodeByClassByName("vtkMRMLMarkupsLineNode", lineNodeName)
        lineNode.AddObserver(slicer.vtkMRMLMarkupsLineNode.PointModifiedEvent, ShowAngle)
    



"""
Created on Mon Apr 17 14:08:38 2023

@author: Dan
Function to project the CSA Angle measured in the 3D slicer model to a 2D angle as in a true-anteroposterior radiograph

"""

# Import CSA Markup


def CSA():

    import math

    ScapulaEnd = slicer.util.getNode("Scapula End")

    GleneDown = slicer.util.getNode("Glene Down")
    GleneUp = slicer.util.getNode("Glene Up")

    # # Get frontal transform node
    # TransformFrontal = slicer.util.getNode("Glene - Rotation plan frontal")

    P1 = np.zeros(3)
    ScapulaEnd.GetNthControlPointPositionWorld(0, P1)

    P2 = np.zeros(3)
    GleneDown.GetNthControlPointPositionWorld(0, P2)

    P3 = np.zeros(3)
    GleneUp.GetNthControlPointPositionWorld(0, P3)

    P1[2], P2[2], P3[2] = 0, 0, 0

    # Vecteurs directeurs
    V1 = P1 - P2
    V2 = P3 - P2

    Angle = 0
    Angle = math.acos(np.dot(V1.T, V2) / (np.linalg.norm(V1) * np.linalg.norm(V2)))
    Angle = Angle / math.pi * 180

    print("CSA = {0:0.3f}".format(Angle) + "°")


"""
CSA variation script
Rotates the glenoid implant around an axis and calculates the matrix rotation required : Variation CSA - Matrice de rotation glene
Outputs the inclination angle and the CSA


In 3d slicer : Modify the IS rotation of the transform : Variation CSA - Rotation glene
"""

GlenRotationAxisName = "Axe anteroposterieur glenoid implant"

# Nom de l'axe de rotation autour duquel tourne l'implant glénoïdien
rotationAxisMarkupsNode = slicer.util.getNode(GlenRotationAxisName)

# # Get frontal transform node
# TransformFrontal = slicer.util.getNode("Glene - Rotation plan frontal")

# This transform can be edited in Transforms module (Edit / Rotation / IS slider)
rotationTransformNode = slicer.util.getNode("Variation CSA - Rotation glene")

# This transform can be edited in Transforms module (Edit / Rotation / LR slider)
AcromionOffsetTransformNode = slicer.util.getNode("Variation CSA - Acromion Offset")

def CSAVariation(unusedArg1=None, unusedArg2=None, unusedArg3=None):
    finalTransformNode = slicer.util.getNode("Variation CSA - Matrice de rotation glene")
    
    
    # Code to rotate a node around an axis (from the Code repositories)
    rotationAxisPoint1_World = np.zeros(3)
    rotationAxisMarkupsNode.GetNthControlPointPositionWorld(0, rotationAxisPoint1_World)
    rotationAxisPoint2_World = np.zeros(3)
    rotationAxisMarkupsNode.GetNthControlPointPositionWorld(1, rotationAxisPoint2_World)

    axisDirectionZ_World = rotationAxisPoint2_World - rotationAxisPoint1_World
    axisDirectionZ_World = axisDirectionZ_World / np.linalg.norm(axisDirectionZ_World)
    
    
    # Get transformation between world coordinate system and rotation axis aligned coordinate system
    worldToRotationAxisTransform = vtk.vtkMatrix4x4()
    
    
    p = vtk.vtkPlaneSource()
    p.SetNormal(axisDirectionZ_World)
    axisOrigin = np.array(p.GetOrigin())
    axisDirectionX_World = np.array(p.GetPoint1()) - axisOrigin
    axisDirectionY_World = np.array(p.GetPoint2()) - axisOrigin
    rotationAxisToWorldTransform = np.row_stack((np.column_stack((axisDirectionX_World, axisDirectionY_World, axisDirectionZ_World, rotationAxisPoint1_World)), (0, 0, 0, 1)))
    rotationAxisToWorldTransformMatrix = slicer.util.vtkMatrixFromArray(rotationAxisToWorldTransform)
    worldToRotationAxisTransformMatrix = slicer.util.vtkMatrixFromArray(np.linalg.inv(rotationAxisToWorldTransform))
    # Compute transformation chain
    rotationMatrix = vtk.vtkMatrix4x4()
    rotationTransformNode.GetMatrixTransformToParent(rotationMatrix)
    finalTransform = vtk.vtkTransform()
    finalTransform.Concatenate(rotationAxisToWorldTransformMatrix)
    finalTransform.Concatenate(rotationMatrix)
    finalTransform.Concatenate(worldToRotationAxisTransformMatrix)
    finalTransformNode.SetMatrixTransformToParent(finalTransform.GetMatrix())

    # shows the inclination after rotation
    ShowAngle("Scapula - Axe transverse", "Glene - Up/Down Axis", "Inclination")

    # Show version
    ShowAngle("Scapula - Axe transverse", "Glene - Left/Right Axis", "Version")

    # Shows the acromion offset
    acromion_offset_matrix = vtk.vtkMatrix4x4()

    # Gets the matrix transform of the Acromion Offset transform
    AcromionOffsetTransformNode.GetMatrixTransformToParent(acromion_offset_matrix)

    # Gets the mediolateral offset
    acromion_offset = acromion_offset_matrix.GetElement(0, 3)

    print("Acromion Lateral Offset = {0:0.3f} mm".format(acromion_offset))
        
    # Shows CSA
    CSA()
    print("\n")
    


# Manual initial update
CSAVariation()

# Automatic update when point is moved or transform is modified
rotationTransformNodeObserver = rotationTransformNode.AddObserver(slicer.vtkMRMLTransformNode.TransformModifiedEvent, CSAVariation)
rotationAxisMarkupsNodeObserver = rotationAxisMarkupsNode.AddObserver(slicer.vtkMRMLMarkupsNode.PointModifiedEvent, CSAVariation)
AcromionOffsetTransformNodeObserver = AcromionOffsetTransformNode.AddObserver(slicer.vtkMRMLMarkupsNode.TransformModifiedEvent, CSAVariation)

# Execute these lines to stop automatic updates:
# rotationTransformNode.RemoveObserver(rotationTransformNodeObserver)
# rotationAxisMarkupsNode.RemoveObserver(rotationAxisMarkupsNodeObserver)
# AcromionOffsetTransformNode.RemoveObserver(AcromionOffsetTransformNodeObserver)

"""
Created on Tue May 23 09:12:43 2023

@author: userscale fator
Function that calculates the CSA based on the acromion Length offset in the x scapula coordinate in mm
"""


def AcromionOffset2CSA(AcromionOffset: float):

    """
    Acromion offset in millimeter
    """
    # Get frontal transform node
    # TransformFrontal = slicer.util.getNode("Glene - Rotation plan frontal")

    GleneDown = slicer.util.getNode("Glene Down")
    GleneUp = slicer.util.getNode("Glene Up")
    ScapulaEnd = slicer.util.getNode("Scapula End")

    P0 = np.zeros(3)
    ScapulaEnd.GetNthControlPointPositionWorld(0, P0)

    P2 = np.zeros(3)
    GleneDown.GetNthControlPointPositionWorld(0, P2)

    P3 = np.zeros(3)
    GleneUp.GetNthControlPointPositionWorld(0, P3)

    Offset = np.array([[AcromionOffset], [0], [0]])
    P1 = P0 + Offset

    # CSA calculation
    P1[2], P2[2], P3[2] = 0, 0, 0

    V1 = P1 - P2
    V2 = P3 - P2

    Angle = math.acos(np.dot(V1.T, V2) / (np.linalg.norm(V1) * np.linalg.norm(V2)))
    Angle = Angle / math.pi * 180

    print("CSA = {0}".format(Angle))


def ExportMat2Anybody(NodeName: str):

    MatrixNode = slicer.util.getNode(NodeName)
    Matrix = vtk.vtkMatrix4x4()
    MatrixNode.GetMatrixTransformToWorld(Matrix)

    Matrix = slicer.util.arrayFromVTKMatrix(Matrix)

    Matrix = np.around(Matrix, decimals=6)

    Rotation = f'AnyMat33 Rotation = {{\n{{{Matrix[0,0]} , {Matrix[0,1]} , {Matrix[0,2]}}},\n{{{Matrix[1,0]} , {Matrix[1,1]} , {Matrix[1,2]}}},\n{{{Matrix[2,0]} , {Matrix[2,1]} , {Matrix[2,2]}}}\n\n}};'
    Translation = f'AnyVec3 Position = 0.001*{{{Matrix[0,3]} , {Matrix[1,3]} , {Matrix[2,3]}}};'

    print(Translation)
    print(Rotation)


"""
Created on Mon May 22 14:32:37 2023

@author: dan
Exports the Glene position and rotation matrix and ghProth Position (center of the sphere fitted to the glenoid implant surface)
"""

def ExportCSA():

    # FrontalTransformOn(True)

    MatrixNode = slicer.util.getNode("GleneImplant - RotMat 4x4")
    ghNode = slicer.util.getNode("ghProth")

    Matrix = vtk.vtkMatrix4x4()
    MatrixNode.GetMatrixTransformToWorld(Matrix)

    Matrix = slicer.util.arrayFromVTKMatrix(Matrix)

    Matrix = np.around(Matrix, decimals=6)

    ghPosition = np.zeros(3)
    ghNode.GetNthControlPointPositionWorld(0, ghPosition)

    ghPosition = np.round(ghPosition, 6)

    inclination_angle = round(get_angle_between_axes("Scapula - Axe transverse", "Glene - Up/Down Axis", "Inclination") - 90, 3)
    version_angle = round(get_angle_between_axes("Scapula - Axe transverse", "Glene - Left/Right Axis", "Inclination") - 90, 3)

    RotationAxisName = f'AnyString RotationAxis = "{GlenRotationAxisName}";'
    Inclination = f"AnyVar GleneImplantTiltAngle = {inclination_angle};"
    Version = f"AnyVar GleneImplantVersionAngle = {version_angle};"
    Rotation = f'AnyMat33 Rotation = {{\n{{{Matrix[0,0]} , {Matrix[0,1]} , {Matrix[0,2]}}},\n{{{Matrix[1,0]} , {Matrix[1,1]} , {Matrix[1,2]}}},\n{{{Matrix[2,0]} , {Matrix[2,1]} , {Matrix[2,2]}}}\n\n}};'
    Translation = f'AnyVec3 Position = 0.001*{{{Matrix[0,3]} , {Matrix[1,3]} , {Matrix[2,3]}}};'
    ghProthLocal = f'\nAnyVec3 Center_Absolute = 0.001*{{{ghPosition[0]} , {ghPosition[1]} , {ghPosition[2]}}};'

    print(RotationAxisName)
    print(Inclination)
    print(Version)
    print(Translation)
    print(Rotation)
    print(ghProthLocal)

# %% Tools for converting Anybody to python
""""
Converts a string containing the anybody definition of a matrix and converts it to an numpy.array
"""


def AnyMatrix2Array(string):
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


"""
Input : rotation matrix and translation vector to a 4x4 matrix for a 3d slicer transform
"""


def TransformMatrix(RotMat=None, Vect=None):

    if RotMat is None:
        RotMat = np.eye(3)
    if Vect is None:
        Vect = np.zeros(3)

    Mat = np.eye(4)

    Mat[0:3, 0:3] = RotMat
    Mat[0:3, 3] = Vect
    return Mat


"""
Function to set the 4x4 transform matrix to a transform in 3d slicer from a rotation matrix and translation vector numpy arrays
"""


def SetTransformMatrix(TransformName, RotMat=None, Vect=None):

    TransformMat = TransformMatrix(RotMat, Vect)

    Transform = getNode(TransformName)
    Transform.SetMatrixTransformToParent(slicer.util.vtkMatrixFromArray(TransformMat))
