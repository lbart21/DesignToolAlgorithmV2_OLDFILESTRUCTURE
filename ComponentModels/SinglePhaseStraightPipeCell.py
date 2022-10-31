import math as m
import importlib
import sys

class SinglePhaseStraightPipeCell():
    def __init__(self, cell_ID, fluidPair, label) -> None:
        self.GEO = {}
        self.flowState = {}
        self.cell_ID = cell_ID
        self.fluidPair = fluidPair
        self.WestInterface = None
        self.EastInterface = None
        self.label = label
        self.phase = "Single"
        self.InteriorCellFlag = True

        self.conservedProperties = {
            "mass"      : 0.0,
            "xMom"      : 0.0,
            "energy"    : 0.0
        }

        [fluidDataFileName, fluidEOSFileName] = fluidPair["fluid"]

        sys.path.insert(0, "/home/s4393747/dgd/examples/LukesExamples/Algorithms/DesignToolAlgorithmV1/FluidModel")
        fluidDataFile = importlib.import_module(fluidDataFileName)
        fluidEOSFile = importlib.import_module(fluidEOSFileName)

        self.flowState = fluidEOSFile.FluidModel(fluidDefinition = fluidDataFile.constants) # Make the flowState an 
                                                                                            # object with constants already defined in it (gamma, Cp, etc), 
                                                                                            # then update is inherent from flowState object, so don't have to reimport every time.
                                                                                            # Same occurs in left and right states of interface.
        
    def fillGeometry(self, Geometry):
        self.GEO = Geometry
        
    def fillProps(self, prop1, prop2, propsGiven, vel_x):
        if propsGiven == "pT":
            self.flowState.fs["p"] = prop1
            self.flowState.fs["T"] = prop2
            self.flowState.fs["vel_x"] = vel_x
            self.flowState.UpdateFromPT() 
        if propsGiven == "rhoP":
            self.flowState.fs["p"] = prop1
            self.flowState.fs["rho"] = prop2
            self.flowState.fs["vel_x"] = vel_x
            self.flowState.UpdateFromRhoP()
        self.initialiseConservedProperties()

    def updatePrimativeProperties(self):
        new_density = self.conservedProperties["mass"] 
        new_vel_x = self.conservedProperties["xMom"] / self.conservedProperties["mass"] 
        new_u = self.conservedProperties["energy"] / self.conservedProperties["mass"] - 0.5 * new_vel_x ** 2.0
        self.flowState.fs["vel_x"] = new_vel_x
        self.flowState.fs["rho"] = new_density
        self.flowState.fs["u"] = new_u
        self.flowState.UpdateFromRhoU()
        
    def initialiseConservedProperties(self):
        self.conservedProperties["mass"] = self.flowState.fs["rho"]
        self.conservedProperties["xMom"] = self.flowState.fs["rho"] * self.flowState.fs["vel_x"]
        self.conservedProperties["energy"] = self.flowState.fs["rho"] * (self.flowState.fs["u"] + 0.5 * self.flowState.fs["vel_x"] ** 2)

    def maxAllowableDt(self, cfl):
        return cfl * self.GEO["dx"] / (abs(self.flowState.fs["vel_x"]) + self.flowState.fs["a"])
    
    def completeCellMethods(self, cellArray):
        self.updatePrimativeProperties()