from copy import deepcopy

def FromStagnationInFlow_BC(mesh, BC):
    """
    [p_stag, T_stag] = BC[1][2]
    Take p_stag and T_stag to create a stagnation state. Use velocity on inside of 
    boundary and stagnation enthalpy to find static enthalpy. Update properties from 
    stagnation entropy and static enthalpy. All cells have the interior cell's velocity.
    Velocity is set to 0 if flow is out of boundary.
    Then update the flow properties to complete the state.
    """
    vel_x_boundary = mesh.cellArray[0].flowState.fs["vel_x"]
    [p_stag, T_stag] = BC[1][2]
    bulk_speed = max(0.0, vel_x_boundary) # 0.0 if vel_x_boundary negative, vel_x_boundary if positive, 
                                          # this assumes flow in is always on
                                          # the west boundary
    flowState = deepcopy(mesh.cellArray[0].flowState)
    flowState.fs["p"] = p_stag
    flowState.fs["T"] = T_stag
    flowState.fs["vel_x"] = bulk_speed
    flowState.UpdateFromPT()
    stagnation_enthalpy = flowState.fs["h"]
    static_enthalpy = stagnation_enthalpy - 0.5 * bulk_speed ** 2.0
    flowState.fs["h"] = static_enthalpy
    flowState.UpdateFromHS()

    for cell in range(len(mesh.cellArray)):
        mesh.cellArray[cell].flowState = flowState
        
    return mesh