from copy import deepcopy

def FromStagnationWithMassFlowRateInFlow_BC(mesh, BC):
    """
    [p_stag, T_stag] = BC[1][2]
    Take p_stag and T_stag to create a stagnation state. Use velocity on inside of 
    boundary and stagnation enthalpy to find static enthalpy. Update properties from 
    stagnation entropy and static enthalpy. All cells have the interior cell's velocity.
    Velocity is set to 0 if flow is out of boundary.
    Then update the flow properties to complete the state.
    """
    relaxation_factor = 0.1
    [p_stag, T_stag, mass_flux] = BC[1][2]
    p0_min = 0.1 * p_stag
    p0_max = 10.0 * p_stag
    vel_x_boundary = mesh.cellArray[0].flowState.fs["vel_x"]
    p_boundary = mesh.cellArray[0].flowState.fs["p"]
    rho_boundary = mesh.cellArray[0].flowState.fs["rho"]
    rhoU_boundary = rho_boundary * vel_x_boundary
    A_boundary = mesh.interfaceArray[0].GEO["A"]
    dp_over_p = 0.5 * relaxation_factor / rho_boundary * ( (mass_flux / A_boundary) ** 2.0 \
                                                            - rhoU_boundary * abs(rhoU_boundary) ) / p_boundary
    new_p0 = (1.0 + dp_over_p) * p_stag
    new_p0 = min(max(new_p0, p0_min), p0_max)
    flowState = deepcopy(mesh.cellArray[0].flowState)
    flowState.fs["p"] = new_p0
    flowState.fs["T"] = T_stag
    flowState.UpdateFromPT()
    stagnation_enthalpy = flowState.fs["h"]
    bulk_speed = max(0.0, vel_x_boundary) # 0.0 if vel_x_boundary negative, vel_x_boundary if positive, 
                                          # this assumes flow in is always on
                                          # the west boundary
    static_enthalpy = stagnation_enthalpy - 0.5 * bulk_speed ** 2.0
    
    flowState.fs["vel_x"] = bulk_speed
    
    flowState.fs["h"] = static_enthalpy
    flowState.UpdateFromHS()

    for cell in range(len(mesh.cellArray)):
        mesh.cellArray[cell].flowState = flowState
        
    return mesh