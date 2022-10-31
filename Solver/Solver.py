from DesignToolAlgorithmV1.Solver.WriteToDataFile import WriteToDataFile
from DesignToolAlgorithmV1.Integrate.Integrate import Integrate
import objsize
from copy import deepcopy
class Solver():
    def __init__(self, meshObject, cfl_flag, tFinal, dataSaveDt) -> None:
        """
        meshObject = object with attributes: cellArray, interfaceArray, mapCellIDToWestInterfaceIdx, 
        cfl_flag = [Bool, float]
        labels = 
        """
        tCurrent = 0.0
        self.totalDataDict = {} #We'll store all data for performance calculations but will only write to text files at certain times
        #self.addToData(time = tCurrent, data = meshObject.cellArray)
        WriteToDataFile(cellArray = meshObject.cellArray, time = tCurrent, labels = meshObject.componentLabels)
        
        time_tol = 1e-9
        tWriteTol = 1e-9
        tWrite = dataSaveDt
        writtenData = False
        currentMeshObject = deepcopy(meshObject)
        currentStep = 0
        while tCurrent < tFinal and abs(tCurrent - tFinal) > time_tol:
            newData = Integrate(mesh = currentMeshObject, cfl_flag = cfl_flag, tCurrent = tCurrent, currentStep = currentStep)
            tCurrent += newData.dtTotal
            #print("t: ", tCurrent)
            writtenData = False
            if tCurrent > tWrite - tWriteTol:
                print("Writing data, t = ", tCurrent)
                WriteToDataFile(cellArray = newData.mesh.cellArray, time = tCurrent, labels = newData.mesh.componentLabels)
                writtenData = True
                tWrite += dataSaveDt
            #self.addToData(time = tCurrent, data = newData.mesh.cellArray)
            currentMeshObject = deepcopy(newData.mesh)
            currentStep += 1
            
        if not writtenData:
            WriteToDataFile(cellArray = currentMeshObject.cellArray, time = tCurrent, labels = currentMeshObject.componentLabels)
        
    def addToData(self, time, data):
        self.totalDataDict[str(time)] = data
        print("totalDataDict size:", objsize.get_deep_size(self.totalDataDict) / 1e6, "Mb")