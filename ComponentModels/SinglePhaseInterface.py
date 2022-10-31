import math as m
import importlib
import sys
from DesignToolAlgorithmV1.Reconstruction.Reconstruction import getReconstruction
from DesignToolAlgorithmV1.Fluxes.fluidFluxes import *
from DesignToolAlgorithmV1.Reconstruction.LocateNeighbouringCellIndices import findIdxOfCellRecursively

class SinglePhaseInterface():
    def __init__(self, interface_ID, nL, nR, fluidPair, reconstructionScheme, limiter, reconstructionProperties, updateFrom, fluxScheme) -> None:
        self.interface_ID = interface_ID
        self.nL = nL
        self.nR = nR
        self.fluidPair = fluidPair
        self.fluxScheme = fluxScheme
        self.FluxFlag = True
        self.reconstructionScheme = reconstructionScheme
        self.limiter = limiter
        self.reconstructionProperties = reconstructionProperties
        self.updateFrom = updateFrom
        self.stencilBeenFormed = False
        
        [fluidDataFileName, fluidEOSFileName] = fluidPair["fluid"]
        sys.path.insert(0, "/home/s4393747/dgd/examples/LukesExamples/Algorithms/DesignToolAlgorithmV1/FluidModel")
        fluidDataFile = importlib.import_module(fluidDataFileName)
        fluidEOSFile = importlib.import_module(fluidEOSFileName)

        self.LftState = fluidEOSFile.FluidModel(fluidDefinition = fluidDataFile.constants)
        self.RghtState = fluidEOSFile.FluidModel(fluidDefinition = fluidDataFile.constants)

    def completeInterfaceMethods(self, cellArray, interfaceArray, mapCellIDToWestInterfaceIdx, \
                                        mapCellIDToEastInterfaceIdx, mapInterfaceIDToWestCellIdx, \
                                        mapInterfaceIDToEastCellIdx, dt_inv):
        if self.FluxFlag:
            if not self.stencilBeenFormed: 
                self.FormListsOfStencilIndices(mapCellIDToWestInterfaceIdx = mapCellIDToWestInterfaceIdx, \
                                                mapCellIDToEastInterfaceIdx = mapCellIDToEastInterfaceIdx, \
                                                mapInterfaceIDToWestCellIdx = mapInterfaceIDToWestCellIdx, \
                                                mapInterfaceIDToEastCellIdx = mapInterfaceIDToEastCellIdx, \
                                                cellArray = cellArray, interfaceArray = interfaceArray)
                self.stencilBeenFormed = True

            self.reconstructStates(cellArray = cellArray)
            self.calculateFluxes()
            self.updateNeighbouringCellsCPs(dt_inv = dt_inv, cellArray = cellArray)
    
    def fillGeometry(self, Geometry):
        self.GEO = Geometry

    def reconstructStates(self, cellArray):
        """
        reconstructionScheme == ["Copy", [rnL, rnR]] etc etc where rnL and rnR are the REQUIRED stencil lengths in the left and right directions for the chosen scheme
        reconstructionProperties = list of properties to reconstruct eg ["p", "T", "vel_x", "alpha_g"]
        updateFrom = str(fluidPropertyUpdatePair) eg "pT", "rhoP", "rhoU"
        """
        #print(self.interface_ID, self.LftStencilIdxs, self.RghtStencilIdxs)
        if self.FluxFlag == True: #Do reconstruction
            for prop in self.reconstructionProperties:
                qL_stencil, dxL_stencil, \
                qR_stencil, dxR_stencil = self.getStencils(property = prop, cellArray = cellArray)
                LftProp, RghtProp = getReconstruction(reconstruction = self.reconstructionScheme[0], limiter = self.limiter, \
                                                        qL_stencil = qL_stencil, dxL_stencil = dxL_stencil, \
                                                        qR_stencil = qR_stencil, dxR_stencil = dxR_stencil)
                self.LftState.fs[prop] = LftProp
                self.RghtState.fs[prop] = RghtProp
            #print(self.LftState, self.RghtState)
            if self.updateFrom == "pT":
                self.LftState.UpdateFromPT()
                self.RghtState.UpdateFromPT()
            #print(self.LftState.fs["p"], self.RghtState.fs["p"])
        else: #Not an interface to calculate fluxes over, so skip
            pass

    def calculateFluxes(self):
        if self.FluxFlag == True:
            if self.fluxScheme == "AUSMPlusUPOriginal":
                self.boundaryFluxes = AUSMPlusupORIGINAL(   a_L_forMa = self.LftState.fs["a"],      a_R_forMa = self.RghtState.fs["a"], \
                                                            p_L = self.LftState.fs["p"],            rho_L = self.LftState.fs["rho"], \
                                                            h_L = self.LftState.fs["h"],            a_L = self.LftState.fs["a"], \
                                                            vel_x_L = self.LftState.fs["vel_x"],    p_R = self.RghtState.fs["p"], \
                                                            rho_R = self.RghtState.fs["rho"],       h_R = self.RghtState.fs["h"], \
                                                            a_R = self.RghtState.fs["a"],           vel_x_R = self.RghtState.fs["vel_x"]).fluxes

            elif self.fluxScheme == "AUSMPlusUPPaper":
                self.boundaryFluxes = AUSMPlusupPAPER(  a_L_forMa = self.LftState.fs["a"],  a_L = self.LftState.fs["a"], \
                                                        p_L = self.LftState.fs["p"],        h_L = self.LftState.fs["h"], \
                                                        rho_L = self.LftState.fs["rho"],    vel_x_L = self.LftState.fs["vel_x"], \
                                                        a_R_forMa = self.RghtState.fs["a"], a_R = self.RghtState.fs["a"], \
                                                        p_R = self.RghtState.fs["p"],       h_R = self.RghtState.fs["h"], \
                                                        rho_R = self.RghtState.fs["rho"],   vel_x_R = self.RghtState.fs["vel_x"]).fluxes
            
        else:
            pass

    def updateNeighbouringCellsCPs(self, dt_inv, cellArray):
        if self.FluxFlag == True:
            westCell = cellArray[self.LftStencilIdxs[0]]
            if westCell.InteriorCellFlag: # Check if cell is a ghost cell as we don't want to update the conserved properties of ghost cells
                westCell.conservedProperties["mass"] -= dt_inv * self.GEO["A"] * self.boundaryFluxes["mass"] / westCell.GEO["dV"]
                westCell.conservedProperties["xMom"] -= dt_inv * self.GEO["A"] * self.boundaryFluxes["xMom"] / westCell.GEO["dV"]
                westCell.conservedProperties["xMom"] -= dt_inv * westCell.GEO["A_c"] * self.boundaryFluxes["p"] / westCell.GEO["dV"]
                westCell.conservedProperties["energy"] -= dt_inv * self.GEO["A"] * self.boundaryFluxes["energy"] / westCell.GEO["dV"]

            eastCell = cellArray[self.RghtStencilIdxs[0]]
            if eastCell.InteriorCellFlag: # Check if cell is a ghost cell as we don't want to update the conserved properties of ghost cells
                eastCell.conservedProperties["mass"] += dt_inv * self.GEO["A"] * self.boundaryFluxes["mass"] / eastCell.GEO["dV"]
                eastCell.conservedProperties["xMom"] += dt_inv * self.GEO["A"] * self.boundaryFluxes["xMom"] / eastCell.GEO["dV"]
                eastCell.conservedProperties["xMom"] += dt_inv * eastCell.GEO["A_c"] * self.boundaryFluxes["p"] / eastCell.GEO["dV"]
                eastCell.conservedProperties["energy"] += dt_inv * self.GEO["A"] * self.boundaryFluxes["energy"] / eastCell.GEO["dV"]
                
        else:
            pass

    def getStencils(self, property, cellArray):
        ### Fill left stencil
        qL_stencil = [None] * len(self.LftStencilIdxs)
        dxL_stencil = [None] * len(self.LftStencilIdxs)
        for ind, LftCellIdx in enumerate(self.LftStencilIdxs):
            cell = cellArray[LftCellIdx]
            qL_stencil[ind] = cell.flowState.fs[property]
            dxL_stencil[ind] = cell.GEO["dx"]

        ### Fill right stencil
        qR_stencil = [None] * len(self.RghtStencilIdxs)
        dxR_stencil = [None] * len(self.RghtStencilIdxs)
        for ind, RghtCellIdx in enumerate(self.RghtStencilIdxs):
            cell = cellArray[RghtCellIdx]
            qR_stencil[ind] = cell.flowState.fs[property]
            dxR_stencil[ind] = cell.GEO["dx"]
            
        return qL_stencil, dxL_stencil, qR_stencil, dxR_stencil
        

    def FormListsOfStencilIndices(self, mapCellIDToWestInterfaceIdx, mapCellIDToEastInterfaceIdx, \
                                        mapInterfaceIDToWestCellIdx, mapInterfaceIDToEastCellIdx, cellArray, interfaceArray):
        self.LftStencilIdxs = [None] * self.reconstructionScheme[1][0]
        self.RghtStencilIdxs = [None] * self.reconstructionScheme[1][1]
        for LftStencilIdx in range(self.reconstructionScheme[1][0]):
            cellIdx = findIdxOfCellRecursively(interfaceID = self.interface_ID, recursionDepth = LftStencilIdx + 1, \
                                                direction = "West", mapInterfaceIDToEastCellIdx = mapInterfaceIDToEastCellIdx, \
                                                cellArray = cellArray, mapCellIDToEastInterfaceIdx = mapCellIDToEastInterfaceIdx, \
                                                mapInterfaceIDToWestCellIdx = mapInterfaceIDToWestCellIdx, \
                                                mapCellIDToWestInterfaceIdx = mapCellIDToWestInterfaceIdx, \
                                                interfaceArray = interfaceArray)
            self.LftStencilIdxs[LftStencilIdx] = cellIdx #In order of [L_Idx0, L_Idx1, ...]
        
        for RghtStencilIdx in range(self.reconstructionScheme[1][1]):
            cellIdx = findIdxOfCellRecursively(interfaceID = self.interface_ID, recursionDepth = RghtStencilIdx + 1, \
                                                direction = "East", mapInterfaceIDToEastCellIdx = mapInterfaceIDToEastCellIdx, \
                                                cellArray = cellArray, mapCellIDToEastInterfaceIdx = mapCellIDToEastInterfaceIdx, \
                                                mapInterfaceIDToWestCellIdx = mapInterfaceIDToWestCellIdx, \
                                                mapCellIDToWestInterfaceIdx = mapCellIDToWestInterfaceIdx, \
                                                interfaceArray = interfaceArray)
            self.RghtStencilIdxs[RghtStencilIdx] = cellIdx #In order of [R_Idx0, R_Idx1, ...]
   