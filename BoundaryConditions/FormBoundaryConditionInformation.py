from DesignToolAlgorithmV1.ComponentModels.basicMeshObject import basic1DMeshObject
from DesignToolAlgorithmV1.Reconstruction.LocateNeighbouringCellIndices import findIdxOfCellRecursively
from DesignToolAlgorithmV1.Reconstruction.LocateNeighbouringInterfaceIndices import findIdxOfInterfaceRecursively
from DesignToolAlgorithmV1.BoundaryConditions.SimpleOutFlow_BC import SimpleOutFlow_BC
from DesignToolAlgorithmV1.BoundaryConditions.SupersonicInflow_BC import SupersonicInFlow_BC
from DesignToolAlgorithmV1.BoundaryConditions.FromStagnationInFlow_BC import FromStagnationInFlow_BC
from DesignToolAlgorithmV1.BoundaryConditions.FromStagnationWithMassFlowRateInFlow_BC import FromStagnationWithMassFlowRateInFlow_BC
from DesignToolAlgorithmV1.BoundaryConditions.FixedPOutFlow_BC import FixedPOutFlow_BC
from DesignToolAlgorithmV1.BoundaryConditions.FixedPTOutFlow_BC import FixedPTOutFlow_BC
from DesignToolAlgorithmV1.BoundaryConditions.WallWithSlip_BC import WallWithSlip_BC
from DesignToolAlgorithmV1.BoundaryConditions.WallNoSlip_BC import WallNoSlip_BC
from copy import deepcopy

class FormBoundaryConditionInformation():
    def __init__(self, mesh, BC) -> None:
        ### Find out if BC is on east or west boundary
        onWestBoundaryBool = False
        onEastBoundaryBool = False
        if mesh.mapInterfaceIDToWestCellIdx[BC[0]] == None:
            onWestBoundaryBool = True

        elif mesh.mapInterfaceIDToEastCellIdx[BC[0]] == None:
            onEastBoundaryBool = True
        
        else:
            print("Boundary is not on an edge of the mesh")
        ### Wall BCs
        if BC[1][0] == "WallNoSlip_BC":
            if onEastBoundaryBool:
                self.ghostCellLayer = basic1DMeshObject(nCells = BC[1][1], reversed = False)
            elif onWestBoundaryBool:
                self.ghostCellLayer = basic1DMeshObject(nCells = BC[1][1], reversed = True)
            # Copy cells from inner layer to get correct geometry props 
            self.MirrorCopyCells(BC = BC, onEastBoundaryBool = onEastBoundaryBool, \
                                    onWestBoundaryBool = onWestBoundaryBool, mesh = mesh)
            self.ghostCellLayer = WallNoSlip_BC(mesh = self.ghostCellLayer)
        
        elif BC[1][0] == "WallWithSlip_BC":
            if onEastBoundaryBool:
                self.ghostCellLayer = basic1DMeshObject(nCells = BC[1][1], reversed = False)
            elif onWestBoundaryBool:
                self.ghostCellLayer = basic1DMeshObject(nCells = BC[1][1], reversed = True)
            # Copy cells from inner layer to get correct geometry props 
            self.MirrorCopyCells(BC = BC, onEastBoundaryBool = onEastBoundaryBool, \
                                    onWestBoundaryBool = onWestBoundaryBool, mesh = mesh)
            self.ghostCellLayer = WallWithSlip_BC(mesh = self.ghostCellLayer)

        ### InFlow BCs        
        elif BC[1][0] == "SupersonicInFlow_BC":
            if onEastBoundaryBool:
                self.ghostCellLayer = basic1DMeshObject(nCells = BC[1][1], reversed = False)
            elif onWestBoundaryBool:
                self.ghostCellLayer = basic1DMeshObject(nCells = BC[1][1], reversed = True)
            # Copy cells from inner layer to get correct geometry props 
            self.MirrorCopyCells(BC = BC, onEastBoundaryBool = onEastBoundaryBool, \
                                    onWestBoundaryBool = onWestBoundaryBool, mesh = mesh)
            self.ghostCellLayer = SupersonicInFlow_BC(mesh = self.ghostCellLayer, BC = BC)
        
        elif BC[1][0] == "FromStagnationInFlow_BC":
            if onEastBoundaryBool:
                self.ghostCellLayer = basic1DMeshObject(nCells = BC[1][1], reversed = False)
            elif onWestBoundaryBool:
                self.ghostCellLayer = basic1DMeshObject(nCells = BC[1][1], reversed = True)
            # Copy cells from inner layer to get correct geometry props 
            self.MirrorCopyCells(BC = BC, onEastBoundaryBool = onEastBoundaryBool, \
                                    onWestBoundaryBool = onWestBoundaryBool, mesh = mesh)
            self.ghostCellLayer = FromStagnationInFlow_BC(mesh = self.ghostCellLayer, BC = BC)
    
        elif BC[1][0] == "FromStagnationWithMassFlowRateInFlow_BC":
            if onEastBoundaryBool:
                self.ghostCellLayer = basic1DMeshObject(nCells = BC[1][1], reversed = False)
            elif onWestBoundaryBool:
                self.ghostCellLayer = basic1DMeshObject(nCells = BC[1][1], reversed = True)
            # Copy cells from inner layer to get correct geometry props 
            self.MirrorCopyCells(BC = BC, onEastBoundaryBool = onEastBoundaryBool, \
                                    onWestBoundaryBool = onWestBoundaryBool, mesh = mesh)
            self.ghostCellLayer = FromStagnationWithMassFlowRateInFlow_BC(mesh = self.ghostCellLayer, BC = BC)
            
        ### OutFlow BCs
        elif BC[1][0] == "SimpleOutFlow_BC":
            if onEastBoundaryBool:
                self.ghostCellLayer = basic1DMeshObject(nCells = BC[1][1], reversed = False)
            elif onWestBoundaryBool:
                self.ghostCellLayer = basic1DMeshObject(nCells = BC[1][1], reversed = True)

            # Copy cells and interfaces from inner layer to get correct geometry props 
            self.MirrorCopyCells(BC = BC, onEastBoundaryBool = onEastBoundaryBool, \
                                    onWestBoundaryBool = onWestBoundaryBool, mesh = mesh)
            self.ghostCellLayer = SimpleOutFlow_BC(mesh = self.ghostCellLayer)

        elif BC[1][0] == "SimpleExtrapolateOutFlow_BC":
            pass

        elif BC[1][0] == "FixedPOutFlow_BC":
            if onEastBoundaryBool:
                self.ghostCellLayer = basic1DMeshObject(nCells = BC[1][1], reversed = False)
            elif onWestBoundaryBool:
                self.ghostCellLayer = basic1DMeshObject(nCells = BC[1][1], reversed = True)

            # Copy cells and interfaces from inner layer to get correct geometry props 
            self.MirrorCopyCells(BC = BC, onEastBoundaryBool = onEastBoundaryBool, \
                                    onWestBoundaryBool = onWestBoundaryBool, mesh = mesh)
            self.ghostCellLayer = FixedPOutFlow_BC(mesh = self.ghostCellLayer, BC = BC)

        elif BC[1][0] == "FixedPTOutFlow_BC":
            if onEastBoundaryBool:
                self.ghostCellLayer = basic1DMeshObject(nCells = BC[1][1], reversed = False)
            elif onWestBoundaryBool:
                self.ghostCellLayer = basic1DMeshObject(nCells = BC[1][1], reversed = True)

            # Copy cells and interfaces from inner layer to get correct geometry props 
            self.MirrorCopyCells(BC = BC, onEastBoundaryBool = onEastBoundaryBool, \
                                    onWestBoundaryBool = onWestBoundaryBool, mesh = mesh)
            self.ghostCellLayer = FixedPTOutFlow_BC(mesh = self.ghostCellLayer, BC = BC)
        


    def MirrorCopyCells(self, BC, onEastBoundaryBool, onWestBoundaryBool, mesh):
        for cell in range(BC[1][1]): #Index [1][1] of BC gives how many ghost cells get generated
            if onEastBoundaryBool:
                cellIdx = findIdxOfCellRecursively(interfaceID = BC[0], recursionDepth = cell + 1, \
                                                        direction = "West", mapInterfaceIDToEastCellIdx = mesh.mapInterfaceIDToEastCellIdx, \
                                                        mapCellIDToEastInterfaceIdx = mesh.mapCellIDToEastInterfaceIdx, \
                                                        mapInterfaceIDToWestCellIdx = mesh.mapInterfaceIDToWestCellIdx, \
                                                        mapCellIDToWestInterfaceIdx = mesh.mapCellIDToWestInterfaceIdx, \
                                                        cellArray = mesh.cellArray, interfaceArray = mesh.interfaceArray)
            elif onWestBoundaryBool:
                cellIdx = findIdxOfCellRecursively(interfaceID = BC[0], recursionDepth = cell + 1, \
                                                        direction = "East", mapInterfaceIDToEastCellIdx = mesh.mapInterfaceIDToEastCellIdx, \
                                                        mapCellIDToEastInterfaceIdx = mesh.mapCellIDToEastInterfaceIdx, \
                                                        mapInterfaceIDToWestCellIdx = mesh.mapInterfaceIDToWestCellIdx, \
                                                        mapCellIDToWestInterfaceIdx = mesh.mapCellIDToWestInterfaceIdx, \
                                                        cellArray = mesh.cellArray, interfaceArray = mesh.interfaceArray)
            ghostCell = deepcopy(mesh.cellArray[cellIdx])
            ghostCell.InteriorFlag = False
            self.ghostCellLayer.cellArray[cell] = ghostCell
            
        for interface in range(BC[1][1] + 1):
            if onEastBoundaryBool:
                interfaceIdx = findIdxOfInterfaceRecursively(interfaceID = BC[0], recursionDepth = interface, \
                                                                direction = "West", mapInterfaceIDToEastCellIdx = mesh.mapInterfaceIDToEastCellIdx, \
                                                                mapCellIDToEastInterfaceIdx = mesh.mapCellIDToEastInterfaceIdx, \
                                                                mapInterfaceIDToWestCellIdx = mesh.mapInterfaceIDToWestCellIdx, \
                                                                mapCellIDToWestInterfaceIdx = mesh.mapCellIDToWestInterfaceIdx, \
                                                                cellArray = mesh.cellArray, interfaceArray = mesh.interfaceArray)
            elif onWestBoundaryBool:
                interfaceIdx = findIdxOfInterfaceRecursively(interfaceID = BC[0], recursionDepth = interface, \
                                                                direction = "East", mapInterfaceIDToEastCellIdx = mesh.mapInterfaceIDToEastCellIdx, \
                                                                mapCellIDToEastInterfaceIdx = mesh.mapCellIDToEastInterfaceIdx, \
                                                                mapInterfaceIDToWestCellIdx = mesh.mapInterfaceIDToWestCellIdx, \
                                                                mapCellIDToWestInterfaceIdx = mesh.mapCellIDToWestInterfaceIdx, \
                                                                cellArray = mesh.cellArray, interfaceArray = mesh.interfaceArray)
            ghostInterface = deepcopy(mesh.interfaceArray[interfaceIdx])
            ghostInterface.FluxFlag = False
            self.ghostCellLayer.interfaceArray[interface] = ghostInterface
        
