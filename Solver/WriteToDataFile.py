import os

def WriteToDataFile(cellArray, time, labels) -> None:
    cwd = os.getcwd()
    for compLabel in labels:
        filename = "dataAt" + str(format(time, ".9f")) +"forComponent" + compLabel + ".txt"
        file = open(cwd + "/data/" + filename, "w")
        file.write("Label: " + compLabel + "\n")
        file.write("Time: " + str(time) + "\n")
        nCells = len(cellArray)
        firstCellInComponent = True
        for cell_idx in range(nCells):
            if cellArray[cell_idx].label == compLabel:
                if cellArray[cell_idx].phase == "Single":
                    jointDataForCell = {**cellArray[cell_idx].flowState.fs, **cellArray[cell_idx].GEO}
                    pass
                else: #Two phase cell, add "_g" to all keys in gas flowState and "_l" to keys in liquid flowState
                    for key in cellArray[cell_idx].gasFlowState.fs.keys():
                        cellArray[cell_idx].gasFlowState.fs[key + "_g"] = cellArray[cell_idx].gasFlowState.fs.pop(key)
                    for key in cellArray[cell_idx].liquidFlowState.fs.keys():
                        cellArray[cell_idx].liquidFlowState.fs[key + "_l"] = cellArray[cell_idx].liquidFlowState.fs.pop(key)
                    jointDataForCell = {**cellArray[cell_idx].gasFlowState.fs, **cellArray[cell_idx].liquidFtate.fs, **cellArray[cell_idx].GEO}
                    
                variableNames = list(jointDataForCell.keys())
                if firstCellInComponent == True:
                    # Add variable names to file
                    file.write("Variables: " + str(len(variableNames)) + "\n")
                    file.write(" ".join(variableNames) + "\n")
                    firstCellInComponent = False
                for name in variableNames:
                    file.write(str(format(jointDataForCell[name], ".9f")) + " ")
                file.write("\n")
        file.close()