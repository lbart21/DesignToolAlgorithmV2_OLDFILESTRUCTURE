from DesignToolAlgorithmV1.ComponentModels.SinglePhaseStraightPipeCell import SinglePhaseStraightPipeCell
from DesignToolAlgorithmV1.ComponentModels.SinglePhaseInterface import SinglePhaseInterface
from DesignToolAlgorithmV1.ComponentModels.meshObject import meshObject

import numpy as np

class SinglePhaseStraightPipe():
    def __init__(self, Geometry, prop1Bounds, prop2Bounds, propsGiven, vel_xBounds, fluidFile, fluidModel, interfaceIndex, \
                        reconstructionScheme, limiter, reconstructionProperties, updateFrom, componentLabel, fluxScheme) -> None:
        """
        Geometry = [D, L, nCells]
        """
        nCells = Geometry[2]
        self.meshObject = meshObject(nCells = nCells)
                              # In run file, set boundary condition by [[Int_ID, [BC_defintion]], [Int_ID, [BC_defintion]], ...]
                                                            # When joining components, need to do mapping from original component boundary condition interface ID to new interface ID
                                                            # This will be the same as when adding the boundary ghost cells
        fluidPairs = {}
        fluidPairs["fluid"] = [fluidFile, fluidModel]  #Dictionary of list of strings of file names 
                                                            #fluidFile = file name in FluidModel directory that defines parameters which specify the fluid
                                                            #fluidModel = file name in FluidModel directory that defines the equation of state for the fluid
                                                            #data in fluidFile and fluidModel must be compatible (no missing information that's needed by EOS)
        self.meshObject.componentLabels = [componentLabel]
        ### Define cells
        self.initialiseCells(Geometry = Geometry, prop1Bounds = prop1Bounds, \
                                prop2Bounds = prop2Bounds, propsGiven = propsGiven, \
                                interfaceIndex = interfaceIndex, vel_xBounds = vel_xBounds, \
                                componentLabel = componentLabel, fluidPairs = fluidPairs)

        ### Define interfaces
        self.initialiseInterfaces(Geometry = Geometry, \
                                reconstructionScheme = reconstructionScheme, limiter = limiter, \
                                reconstructionProperties = reconstructionProperties, \
                                updateFrom = updateFrom, fluxScheme = fluxScheme, fluidPairs = fluidPairs)

        ### Form connections
        self.connectCellsToInterfaces()

    def initialiseCells(self, Geometry, prop1Bounds, prop2Bounds, propsGiven, vel_xBounds, interfaceIndex, componentLabel, fluidPairs):
        [D, L, nCells] = Geometry
        self.nCells = nCells
        
        for cell in range(nCells):
            cellObject = SinglePhaseStraightPipeCell(cell_ID = cell, fluidPair = fluidPairs, label = componentLabel)
            GEO = {
                "dx"    :   L / nCells,
                "dV"    :   0.25 * np.pi * D ** 2 * L / nCells,
                "A_c"   :   0.25 * np.pi * D ** 2,
                "A_s"   :   np.pi * D * L / nCells,
                "pos_x" :   (0.5 + cell) * L / nCells
            }
            cellObject.fillGeometry(Geometry = GEO)
            if cell < interfaceIndex:
                cellObject.fillProps(prop1 = prop1Bounds[0], prop2 = prop2Bounds[0], \
                                        propsGiven = propsGiven, vel_x = vel_xBounds[0])
            else:
                cellObject.fillProps(prop1 = prop1Bounds[1], prop2 = prop2Bounds[1], \
                                        propsGiven = propsGiven, vel_x = vel_xBounds[1])
            self.meshObject.cellArray[cell] = cellObject

    def initialiseInterfaces(self, Geometry, reconstructionScheme, limiter, reconstructionProperties, updateFrom, fluxScheme, fluidPairs):
        [D, L, nCells] = Geometry
        for interface in range(self.nCells + 1):
            interfaceObject = SinglePhaseInterface(interface_ID = interface, nL = interface, nR = self.nCells - interface, \
                                                    fluidPair = fluidPairs, fluxScheme = fluxScheme, \
                                                    reconstructionScheme = reconstructionScheme, limiter = limiter, \
                                                    reconstructionProperties = reconstructionProperties, updateFrom = updateFrom)
            GEO = {"A"  : 0.25 * np.pi * D ** 2}
            interfaceObject.fillGeometry(Geometry = GEO)
            self.meshObject.interfaceArray[interface] = interfaceObject
        
        self.meshObject.boundaryInterfaceIDs = [0, self.nCells]

    def connectCellsToInterfaces(self):
        self.meshObject.mapCellIDToWestInterfaceIdx = [None] * self.nCells
        self.meshObject.mapCellIDToEastInterfaceIdx = [None] * self.nCells
        self.meshObject.mapInterfaceIDToWestCellIdx = [None] * (self.nCells + 1)
        self.meshObject.mapInterfaceIDToEastCellIdx = [None] * (self.nCells + 1)

        for array_idx in range(len(self.meshObject.cellArray)):
            self.meshObject.mapCellIDToWestInterfaceIdx[array_idx] = array_idx
            self.meshObject.mapCellIDToEastInterfaceIdx[array_idx] = array_idx + 1
            self.meshObject.mapInterfaceIDToWestCellIdx[array_idx + 1] = array_idx
            self.meshObject.mapInterfaceIDToEastCellIdx[array_idx] = array_idx
          
    def addBoundaryConditions(self, BC):
        self.meshObject.boundaryConditions.append(BC)






 
            


        


