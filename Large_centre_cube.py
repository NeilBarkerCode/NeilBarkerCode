def BCCFEA03(X, Y, Z, R, Rc, File):
    # Corner spheres same size, internal spheres larger
    # if 2*Rc > BCCEquation there is an interaction between centre spheres
    # eg if R=2, Rc needs to be less than 2.3mm unless you want centre sphere interactions

    import numpy as np

    # BCCEquation packing equation
    BCCEquation = (4 / (np.sqrt(3))) * R

    # create text file with information in it
    TextFile = open(str(File) + 'BCCFEA03-'
                    + str(X) + '-' + str(Y) + '-' + str(Z) + '-' + str(R) + '-' + str(Rc) + '.txt', 'w')

    # create SpaceClaim .py file
    SpaceClaimFile = open(str(File) + 'BCCFEA03-'
                          + str(X) + '-' + str(Y) + '-' + str(Z) + '-' + str(R) + '-' + str(Rc) + '-' + '.py', 'w')

    # Changes SpaceClaim code to 3D mode
    SpaceClaimFile.write('mode = InteractionMode.Solid' + '\n' +
                         'result = ViewHelper.SetViewMode(mode, None)' + '\n')

    # create the cube
    SpaceClaimFile.write(
        'result = BlockBody.Create(Point.Create(MM(0), MM(0), MM(0)), Point.Create(MM(' + str((X - 1) * BCCEquation)
        + '), MM(' + str((Y - 1) * BCCEquation) + '), MM(' + str((Z - 1) * BCCEquation) +
        ')), ExtrudeType.ForceAdd, None)' + '\n')

    # this creates the SpaceClaim corner spheres
    x1 = np.array([1, 0, 0])
    y1 = np.array([0, 1, 0])
    z1 = np.array([0, 0, 1])

    for n1 in range(X):
        for n2 in range(Y):
            for n3 in range(Z):
                xyz1 = n1 * x1 * BCCEquation + n2 * y1 * BCCEquation + n3 * z1 * BCCEquation
                SpaceClaimFile.write(
                    'SphereBody.Create(Point.Create(MM(' + str(xyz1[0]) + '), MM(' + str(xyz1[1]) + '), MM(' + str(
                        xyz1[2]) + ')), Point.Create(MM(' + str(xyz1[0] + R) + '), MM(' + str(xyz1[1]) + '), MM(' + str(
                        xyz1[2]) + ')), ExtrudeType.ForceCut,None)' + '\n')

    # This section creates centre spheres for BCC structures.
    offset = [0.5]
    for n4 in range(X - 1):
        for n5 in range(Y - 1):
            for n6 in range(Z - 1):
                xyz1 = (offset + n4 * x1) * BCCEquation + n5 * y1 * BCCEquation + n6 * z1 * BCCEquation
                SpaceClaimFile.write(
                    'SphereBody.Create(Point.Create(MM(' + str(xyz1[0]) + '), MM(' + str(xyz1[1]) + '), MM(' + str(
                        xyz1[2]) + ')), Point.Create(MM(' + str(xyz1[0] + Rc) + '), MM(' + str(
                        xyz1[1]) + '), MM(' + str(
                        xyz1[2]) + ')), ExtrudeType.ForceCut,None)' + '\n')

    # closes the file just written
    SpaceClaimFile.close()

    # write model information to a text file
    TextFile.write('The Function used: BCCFEA03' + '\n' + '\n')
    TextFile.write('Variables used : ' + '\n' + 'Number of spheres in the X direction ' + str(X)
                   + '\n' + 'Number of spheres in the Y direction ' + str(Y) + '\n'
                   + 'Number of spheres in the Z direction ' + str(Z) + '\n' + 'The corner sphere radius '
                   + str(R) + 'mm' + '\n' + 'The centre sphere radius ' + str(Rc) + '\n' + '\n')

    # calculating the cube dimensions (CHECK TO SEE IF THIS IS CORRECT)
    XDim = (X - 1) * BCCEquation
    YDim = (Y - 1) * BCCEquation
    ZDim = (Z - 1) * BCCEquation

    TextFile.write('The dimensions of the model are:' + '\n' + 'X = ' + str(XDim) + 'mm' + '\n' + 'Y = ' +
                   str(YDim) + 'mm' + '\n' + 'Z = ' + str(ZDim) + 'mm' + '\n' + '\n')

    # number of spheres and windows
    NumberOfCentreSpheres = (X - 1) * (Y - 1) * (Z - 1)
    NumberOfCornerSpheres = (X - 1) * (Y - 1) * (Z - 1)
    NumberOfCentreCornerWindows = NumberOfCornerSpheres * 8

    # The total number of centre centre windows changes with geometry as they are shared with connecting units.
    # In addition, the external windows need to be calculated separately as the lenz volume needs to be halfed

    if (2 * Rc) > BCCEquation:
        NumberOfCentreCentreWindows = (X * ((Y - 1) * (Z - 1)) + Y * ((X - 1) * (Z - 1)) + Z * ((X - 1) * (Y - 1)))
        ExternalNumberOfCentreCentreWindows = ((X - 1) * (Y - 1) * 2) + ((Y - 1) * (Z - 1) * 2) + (
                (Z - 1) * (X - 1) * 2)
        InternalNumberOfCentreCentreWindows = NumberOfCentreCentreWindows - ExternalNumberOfCentreCentreWindows
    else:
        NumberOfCentreCentreWindows = 0
        ExternalNumberOfCentreCentreWindows = 0
        InternalNumberOfCentreCentreWindows = 0

    TextFile.write('The total number of centre spheres - ' + str(NumberOfCentreSpheres) + '\n')
    TextFile.write('The total number of corner spheres - ' + str(NumberOfCornerSpheres) + '\n')
    TextFile.write('The total number of centre centre windows - ' + str(NumberOfCentreCentreWindows) + '\n')
    TextFile.write('The total number of internal centre centre windows - '
                   + str(InternalNumberOfCentreCentreWindows) + '\n')
    TextFile.write('The total number of external centre centre windows - '
                   + str(ExternalNumberOfCentreCentreWindows) + '\n')
    TextFile.write('The total number of centre corner windows - ' + str(NumberOfCentreCornerWindows) + '\n')

    # calculation of the lens between centre centre sphere interaction

    # Distance between centre centre spheres
    DCentreCentre = BCCEquation
    DCentreCorner = 2 * R

    # if the centre spheres do not overlap it returns a complex value. This If statemtent returns the Centre
    # Centre window radius as a zero if this is the case
    if (2 * Rc) > BCCEquation:
        CentreCentreWindowRadius = float((1 / (2 * DCentreCentre)) * np.sqrt((-DCentreCentre + Rc - Rc)
                                                                             * (-DCentreCentre - Rc + Rc) * (
                                                                                     -DCentreCentre + Rc + Rc)
                                                                             * (DCentreCentre + Rc + Rc)))
    else:
        CentreCentreWindowRadius = 0

    CentreCornerWindowRadius = float((1 / (2 * DCentreCorner)) * np.sqrt((-DCentreCorner + R - Rc)
                                                                         * (-DCentreCorner - R + Rc) * (
                                                                                 -DCentreCorner + R + Rc) * (
                                                                                 DCentreCorner + R + Rc)))

    # finding the volume of the removed lens common to the centre centre spheres.
    # If the centre centre window radius is zero the lenz volume is zero.
    if CentreCentreWindowRadius == 0:
        CentreCentreLensVolume = 0
    else:
        CentreCentreLensVolume = float((np.pi * (4 * Rc + DCentreCentre) * ((2 * Rc - DCentreCentre) ** 2)) / 12)

    TotalCentreCentreLensVolume = float((CentreCentreLensVolume * InternalNumberOfCentreCentreWindows)
                                        + ((CentreCentreLensVolume * ExternalNumberOfCentreCentreWindows) / 2))

    # finding the volume of the removed lens common to the centre corner spheres
    CentreCornerLensVolume = float(np.pi * ((Rc + R - DCentreCorner) ** 2) * (DCentreCorner ** 2 + 2 * DCentreCorner * R
                                                                              - 3 * R ** 2 + 2 * DCentreCorner * Rc + 6
                                                                              * R * Rc - 3 * Rc ** 2) / (
                                               12 * DCentreCorner))
    TotalCentreCornerLensVolume = CentreCornerLensVolume * NumberOfCentreCornerWindows

    # The total lens volume within the material.
    TotalLensVolume = TotalCentreCentreLensVolume + TotalCentreCornerLensVolume

    # Total volume of the voids within the cube
    CubeVoidVolume = (NumberOfCentreSpheres * (4 * np.pi * Rc ** 3) / 3) + \
                     (NumberOfCornerSpheres * (4 * np.pi * R ** 3) / 3) - TotalLensVolume

    # Volume of cube without any pores in mm3
    TotalCubeVolume = float((X - 1) * (Y - 1) * (Z - 1) * BCCEquation ** 3)

    # Volume of porous material
    VolumeOfPorousMaterial = TotalCubeVolume - CubeVoidVolume

    # Porosity of cube
    Porosity = float(CubeVoidVolume / TotalCubeVolume)

    # Information for the text file

    TextFile.write('The centre centre window radius is ' + str(CentreCentreWindowRadius) + 'mm' + '\n')
    TextFile.write('The centre corner window radius is ' + str(CentreCornerWindowRadius) + 'mm' + '\n')
    TextFile.write('The single centre centre lens volume - ' + str(CentreCentreLensVolume) + 'mm^3' + '\n')
    TextFile.write('The single centre corner lens volume - ' + str(CentreCornerLensVolume) + 'mm^3' + '\n')
    TextFile.write('The total lens volume - ' + str(TotalLensVolume) + 'mm^3' + '\n')
    TextFile.write('Total volume of pores - ' + str(CubeVoidVolume) + 'mm^3' + '\n')
    TextFile.write('The volume of the the solid cube - ' + str(TotalCubeVolume) + 'mm^3' + '\n')
    TextFile.write('The volume of the porous material - ' + str(VolumeOfPorousMaterial) + '\n')
    TextFile.write('The volume fraction (porosity) of the cube - ' + str(Porosity) + '\n')

    TextFile.close()


File = str('D:/05 Python/12 graded porosity/01 python_output/')


BCCFEA03(2, 2, 7, 2, 2.2, File)
