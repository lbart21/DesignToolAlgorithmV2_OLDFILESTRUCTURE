def SupersonicInFlow_BC(mesh, BC):
    """
    Simply put specified flow state which is in index [1][2] of BC list in order
    [vel_x, p, T] = BC[1][2]
    Then update the flow properties to complete the state.
    All cells have the same properties
    """
    for cell in range(len(mesh.cellArray)):
        [vel_x, p, T] = BC[1][2]
        mesh.cellArray[cell].flowState.fs["vel_x"] = vel_x
        mesh.cellArray[cell].flowState.fs["p"] = p
        mesh.cellArray[cell].flowState.fs["T"] = T
        mesh.cellArray[cell].flowState.UpdateFromPT()

    return mesh