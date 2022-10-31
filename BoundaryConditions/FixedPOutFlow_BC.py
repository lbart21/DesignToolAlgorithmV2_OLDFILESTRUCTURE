def FixedPOutFlow_BC(mesh, BC):
    """
    Ghost cells are mirrored values of the interior cells. Specify outlet pressure
    and update from PT.
    Cells will have same pressure but not necessarily the same velocity, temperature etc.
    """
    [p] = BC[1][2]

    for cell in range(len(mesh.cellArray)):
        mesh.cellArray[cell].flowState.fs["p"] = p
        mesh.cellArray[cell].flowState.UpdateFromPT()

    return mesh