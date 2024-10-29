# -*- coding: utf-8 -*-
"""
Script to copy-paste into the 3d slicer model Variation_CSA.mrml to vary the CSA of the glenoid implant
"""

import numpy as np
import math


# %% Nodes for the CSA variation

GlenRotationAxisName = "Axe anteroposterieur glenoid implant"

# Name of the axis of rotation used to rotate the glenoid implant
rotationAxisMarkupsNode = slicer.util.getNode(GlenRotationAxisName)

# Name of the transform node used by the user to rotate the glenoid implant
# This transform can be edited in Transforms module (Edit / Rotation / IS slider)
rotationTransformNode = slicer.util.getNode("Variation CSA - Rotation glene")

# Name of the node to move the lateral end of the acromion
# This transform can be edited in Transforms module (Edit / Rotation / LR slider)
AcromionOffsetTransformNode = slicer.util.getNode("Variation CSA - Acromion Offset")


# %% CSA variation functions

def get_angle_between_axes(axis_1, axis_2):
    """
    Function that calculates the angle between two axis in degrees
    Adapted from the function ShowAngle in the 3d slicer script repository
    https://slicer.readthedocs.io/en/latest/developer_guide/script_repository.html

    axis1 : str : name of the axis
    axis2 : str : name of the axis

    return
    -------------------------------------
    angleDeg : angle between the axis1 and axis2 in degree
    """

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
def show_glenoid_implant_angle(axis_1, axis_2, name_angle=""):
    """
    Prints angle between 2 axes and subtracts 90° to get angles in the right
    And prints the name of the angle

    Adapted from the function ShowAngle in the 3d slicer script repository
    https://slicer.readthedocs.io/en/latest/developer_guide/script_repository.html

    axis1 : str : name of the axis
    axis2 : str : name of the axis
    name_angle : str : name of the angle

    return
    -------------------------------------
    prints the angle with its name

    """

    # get the angle between the axes and subtracts 90° to be in the right quadrant for the inclination and version
    Angle = get_angle_between_axes(axis_1, axis_2) - 90

    print(name_angle + " = {0:0.3f}".format(Angle) + "°")

    AxisNodeNames = [axis_1, axis_2]

    # Observe line node changes and execute this function each time the axis move
    for lineNodeName in AxisNodeNames:
        lineNode = slicer.util.getFirstNodeByClassByName("vtkMRMLMarkupsLineNode", lineNodeName)
        lineNode.AddObserver(slicer.vtkMRMLMarkupsLineNode.PointModifiedEvent, show_glenoid_implant_angle)


def CSA():
    """
    Function that calculates and prints the CSA

    The function uses the position of three points :
        Acromion End : most lateral point of the acromion
        Glene Down : most inferior point of the glenoid implant articular surface
        Glene Up : most superior point of the glenoid implant articular surface

    The CSA is the angle between two lines (projected on the frontal plane of the scapula):
        line 1 : Acromion End - Glene Down
        line 2 : Glene Down - Glene Up
    """

    ScapulaEnd = slicer.util.getNode("Acromion End")

    GleneDown = slicer.util.getNode("Glene Down")
    GleneUp = slicer.util.getNode("Glene Up")

    P1 = np.zeros(3)
    ScapulaEnd.GetNthControlPointPositionWorld(0, P1)

    P2 = np.zeros(3)
    GleneDown.GetNthControlPointPositionWorld(0, P2)

    P3 = np.zeros(3)
    GleneUp.GetNthControlPointPositionWorld(0, P3)

    # Set the coordinates to 0 to project the angle on the frontal plane of the scapula
    P1[2], P2[2], P3[2] = 0, 0, 0

    # Vecteurs directeurs
    V1 = P1 - P2
    V2 = P3 - P2

    Angle = 0
    Angle = math.acos(np.dot(V1.T, V2) / (np.linalg.norm(V1) * np.linalg.norm(V2)))
    Angle = Angle / math.pi * 180

    print("CSA = {0:0.3f}".format(Angle) + "°")


def CSAVariation(unusedArg1=None, unusedArg2=None, unusedArg3=None):
    """
    CSA variation script
    Rotates the glenoid implant around an axis and calculates the matrix rotation to apply to the glenoid implant ("Variation CSA - Matrice de rotation glene")

    In 3d slicer : Modify the IS rotation of the transform : "Variation CSA - Rotation glene"
    In 3d slicer : Modify the LR translation of the transform : "Variation CSA - Acromion Offset"

    Code adapted from the 3d slicer script repository https://slicer.readthedocs.io/en/latest/developer_guide/script_repository.html#rotate-a-node-around-a-specified-line

    return
    -----------------------------
    Varies the CSA
    Prints the CSA, the inclination, the version and the acromion offset each time the CSA changes
    """

    # Rotation matrix that will be applied to the glenoid implant to rotate it
    finalTransformNode = slicer.util.getNode("Variation CSA - Matrice de rotation glene")

    # Code to rotate a node around an axis
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
    show_glenoid_implant_angle("Scapula - Axe transverse", "Glene - Up/Down Axis", "Inclination")

    # Show version
    show_glenoid_implant_angle("Scapula - Axe transverse", "Glene - Left/Right Axis", "Version")

    # Shows the acromion offset
    acromion_offset_matrix = vtk.vtkMatrix4x4()

    # Gets the matrix transform of the Acromion Offset transform
    AcromionOffsetTransformNode.GetMatrixTransformToParent(acromion_offset_matrix)

    # Gets the mediolateral offset
    acromion_offset = acromion_offset_matrix.GetElement(0, 3)

    print("Acromion Lateral Offset = {0:0.3f} mm".format(acromion_offset))

    # Shows the CSA and the orientation of the glenoid implant
    CSA()
    print("\n")


# Manual initial update of the CSA
CSAVariation()

# Automatic update of the CSA when point is moved or transform is modified
rotationTransformNodeObserver = rotationTransformNode.AddObserver(slicer.vtkMRMLTransformNode.TransformModifiedEvent, CSAVariation)
rotationAxisMarkupsNodeObserver = rotationAxisMarkupsNode.AddObserver(slicer.vtkMRMLMarkupsNode.PointModifiedEvent, CSAVariation)
AcromionOffsetTransformNodeObserver = AcromionOffsetTransformNode.AddObserver(slicer.vtkMRMLMarkupsNode.TransformModifiedEvent, CSAVariation)

# Execute these lines to stop automatic updates:
# rotationTransformNode.RemoveObserver(rotationTransformNodeObserver)
# rotationAxisMarkupsNode.RemoveObserver(rotationAxisMarkupsNodeObserver)
# AcromionOffsetTransformNode.RemoveObserver(AcromionOffsetTransformNodeObserver)


def ExportCSA():
    """
    Function to export the glenoid implant rotation to anybody. Prints a code to be pasted into a ImplanPositions.any file

    Print
    ---------------------------
    RotationAxisName : Name of the rotation axis (AnyString)
    Inclination : Value of the inclination angle (deg) (AnyVar)
    Version : Value of the version angle (deg) (AnyVar)
    Translation : Position of the glenoid implant in the scapula coordinate sytem (sRel value) (AnyVec3)
    Rotation : Rotation matrix of the glenoid implant in the scapula coordinate sytem (ARel Value) (AnyMat33)
    ghProthLocal : Position of the center of a sphere fitted to the glenoid surface in the scapula coordinate system (AnyVec3)
    """

    MatrixNode = slicer.util.getNode("GleneImplant - RotMat 4x4")
    ghNode = slicer.util.getNode("ghProth")

    Matrix = vtk.vtkMatrix4x4()
    MatrixNode.GetMatrixTransformToWorld(Matrix)

    Matrix = slicer.util.arrayFromVTKMatrix(Matrix)

    Matrix = np.around(Matrix, decimals=6)

    ghPosition = np.zeros(3)
    ghNode.GetNthControlPointPositionWorld(0, ghPosition)

    ghPosition = np.round(ghPosition, 6)

    # Version and inclination angle
    inclination_angle = round(get_angle_between_axes("Scapula - Axe transverse", "Glene - Up/Down Axis") - 90, 3)
    version_angle = round(get_angle_between_axes("Scapula - Axe transverse", "Glene - Left/Right Axis") - 90, 3)

    RotationAxisName = f'AnyString RotationAxis = "{GlenRotationAxisName}";'
    Inclination = f"AnyVar GleneImplantInclinationAngle = {inclination_angle};"
    Version = f"AnyVar GleneImplantVersionAngle = {version_angle};"
    Translation = f'AnyVec3 Position = 0.001*{{{Matrix[0,3]} , {Matrix[1,3]} , {Matrix[2,3]}}};'
    Rotation = f'AnyMat33 Rotation = {{\n{{{Matrix[0,0]} , {Matrix[0,1]} , {Matrix[0,2]}}},\n{{{Matrix[1,0]} , {Matrix[1,1]} , {Matrix[1,2]}}},\n{{{Matrix[2,0]} , {Matrix[2,1]} , {Matrix[2,2]}}}\n\n}};'
    ghProthLocal = f'\nAnyVec3 Center_Absolute = 0.001*{{{ghPosition[0]} , {ghPosition[1]} , {ghPosition[2]}}};'

    print(RotationAxisName)
    print(Inclination)
    print(Version)
    print(Translation)
    print(Rotation)
    print(ghProthLocal)
