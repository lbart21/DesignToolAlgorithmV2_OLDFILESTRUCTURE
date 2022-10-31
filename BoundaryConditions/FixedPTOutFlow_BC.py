def FixedPTOutFlow_BC(mesh, BC):
    """
    Ghost cells are mirrored values of the interior cells. Specify outlet pressure
    and temperature and update from PT.
    Cells will have same pressure and temperature (and thus all thermodynamic properties), 
    but not necessarily the same velocity.
    """
    [p, T] = BC[1][2]

    for cell in range(len(mesh.cellArray)):
        mesh.cellArray[cell].flowState.fs["p"] = p
        mesh.cellArray[cell].flowState.fs["T"] = T
        mesh.cellArray[cell].flowState.UpdateFromPT()

    return mesh