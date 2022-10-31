from regex import T
from DesignToolAlgorithmV1.PostProcessing.DataFileToStructuredData import GenerateDataObject
from DesignToolAlgorithmV1.PostProcessing.ProcessEilmerData import ProcessEilmerData
from DesignToolAlgorithmV1.PostProcessing.SIUnitsDictionary import SIUnits
from DesignToolAlgorithmV1.PostProcessing.Symbols import symbols
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.animation import PillowWriter

import os


class GenerateSinglePlots():
    def __init__(self, dataFile, plotVars) -> None:
        self.dataObject = GenerateDataObject(dataFileName = dataFile)
        
        tFinal = self.dataObject.tFinal
        if "Ma" in plotVars:
            self.dataObject.componentData["Ma"] = self.dataObject.componentData["vel_x"] / self.dataObject.componentData["a"]

        if "p_t" in plotVars:
            gamma = self.dataObject.componentData["gamma"]
            p = self.dataObject.componentData["p"]
            Ma = self.dataObject.componentData["Ma"]
            self.dataObject.componentData["p_t"] = p * (1.0 + 0.5 * (gamma - 1.0) * Ma ** 2.0) ** (gamma / (gamma - 1.0))
            
        if "T_t" in plotVars:
            gamma = self.dataObject.componentData["gamma"]
            T = self.dataObject.componentData["T"]
            Ma = self.dataObject.componentData["Ma"]
            self.dataObject.componentData["T_t"] = T * (1.0 + 0.5 * (gamma - 1.0) * Ma ** 2.0)
            
        for var in plotVars:
            fig = plt.figure(figsize=(15, 5))
            formattedTitleTime = '{:.3f}'.format(tFinal / 1e-6)
            formattedFileNameTime = '{:.9f}'.format(tFinal)
            plt.title("Distribution of " + symbols[var] + " at t = " \
                                                    + formattedTitleTime + r'$\mu$' + "s")
            plt.ylabel(symbols[var] + " (" + SIUnits[var] +")", \
                        rotation = "horizontal", ha = "right")
            plt.xlabel("Position (m)")
            plt.scatter(self.dataObject.componentData["pos_x"], self.dataObject.componentData[var], marker = '.')
            filename = var + " distribution at t = " + formattedFileNameTime + ".jpg"
            plt.grid()
            mng = plt.get_current_fig_manager()
            mng.full_screen_toggle()
            currentDir = os.getcwd()
            plt.savefig(currentDir + "/plots/" + filename, bbox_inches="tight")
            plt.close()
        pass

class GenerateWaterfallPlots():
    def __init__(self, dataFiles, plotVars) -> None:
        dataFromFiles = {}
        t_list = []

        for file in dataFiles:
            data = GenerateDataObject(dataFileName = file)
            if "Ma" in plotVars:
                data.componentData["Ma"] = data.componentData["vel_x"] / data.componentData["a"]
            
            if "p_t" in plotVars:
                gamma = data.componentData["gamma"]
                p = data.componentData["p"]
                Ma = data.componentData["Ma"]
                data.componentData["p_t"] = p * (1.0 + 0.5 * (gamma - 1.0) * Ma ** 2.0) ** (gamma / (gamma - 1.0))
            
            if "T_t" in plotVars:
                gamma = data.componentData["gamma"]
                T = data.componentData["T"]
                Ma = data.componentData["Ma"]
                data.componentData["T_t"] = T * (1.0 + 0.5 * (gamma - 1.0) * Ma ** 2.0)

            dataFromFiles[str(data.tFinal)] = data.componentData
            t_list.append(data.tFinal)

        for var in plotVars:
            fig = plt.figure(figsize=(15, 5))
            plt.title("Distribution of " + symbols[var] + " at multiple time values")
            plt.ylabel(symbols[var] + " (" + SIUnits[var] +")", \
                        rotation = "horizontal", ha = "right")
            plt.xlabel("Position (m)")
            for time in t_list:
                formattedTime = '{:.3f}'.format(time / 1e-6)
                plt.scatter(dataFromFiles[str(time)]["pos_x"], dataFromFiles[str(time)][var], \
                                label = "Distribution at t = " + formattedTime + r'$\mu$' + "s", \
                                marker = ".")
            plt.legend()
            plt.grid()
            filename = var + " distribution at multiple times.jpg"
            mng = plt.get_current_fig_manager()
            mng.full_screen_toggle()
            currentDir = os.getcwd()
            plt.savefig(currentDir + "/plots/" + filename, bbox_inches="tight")
            plt.close()

class generateSingleComponentAnimation():
    def __init__(self, dataFiles, slowDownFactor, plotVars) -> None:
        self.data = {}
        self.dataTimes = [0.0]
        for file in dataFiles:
            dataObject = GenerateDataObject(dataFileName = file)
            self.data[str(dataObject.tFinal)] = dataObject.componentData
            self.dataTimes.append(dataObject.tFinal)
        
        timeStepList = []
        for i in range(len(self.dataTimes)):
            timeStepList.append(self.dataTimes[i+1] - self.dataTimes[i])
        
        for var in plotVars:
            self.fig, self.ax = plt.subplots()

            anim = animation.FuncAnimation(self.fig, self.update, frames = self.dataTimes, \
                                            blit = True, repeat = False)
            writerVideo = animation.PillowWriter(fps = 30)
            fileName = "AnimationOf" + var + ".gif"
            currentDir = os.getcwd()
            anim.save(currentDir + "/plots/" + fileName, writer = writerVideo)

    def update(self, i, fargs):
        (var,) = fargs
        x = self.data[str(i)]["pos_x"]
        y = self.data[str(i)][var]
        self.scat = self.ax.scatter(x, y)
        return self.scat, 

        
class GenerateSinglePlotsFromEilmerData():
    def __init__(self, EilmerDataNames, plotVars) -> None:
        EilmerData = ProcessEilmerData(dataFiles = EilmerDataNames)
        t = EilmerData.tFinal
        formattedTitleTime = '{:.3f}'.format(t / 1e-6)
        formattedFileNameTime = '{:.9f}'.format(t)
        if "Ma" in plotVars:
            EilmerData.componentData["Ma"] = EilmerData.componentData["vel_x"] / EilmerData.componentData["a"]

        if "p_t" in plotVars:
            gamma = 1.4
            p = EilmerData.componentData["p"]
            Ma = EilmerData.componentData["Ma"]
            EilmerData.componentData["p_t"] = p * (1.0 + 0.5 * (gamma - 1.0) * Ma ** 2.0) ** (gamma / (gamma - 1.0))
            
        if "T_t" in plotVars:
            gamma = 1.4
            T = EilmerData.componentData["T"]
            Ma = EilmerData.componentData["Ma"]
            EilmerData.componentData["T_t"] = T * (1.0 + 0.5 * (gamma - 1.0) * Ma ** 2.0)
            
        for var in plotVars:
            fig = plt.figure(figsize=(15, 5))
            plt.title("Eilmer Simulation Distribution of " + symbols[var] + " at t = " \
                                                    + formattedTitleTime + r'$\mu$' + "s")
            plt.ylabel(symbols[var] + " (" + SIUnits[var] +")", \
                        rotation = "horizontal", ha = "right")
            plt.xlabel("Position (m)")
            plt.scatter(EilmerData.componentData["pos_x"], EilmerData.componentData[var], marker = '.')
            filename = var + " distribution at t = " + formattedFileNameTime + "WithEilmerSimulation.jpg"
            plt.grid()
            mng = plt.get_current_fig_manager()
            mng.full_screen_toggle()
            currentDir = os.getcwd()
            plt.savefig(currentDir + "/plots/" + filename, bbox_inches="tight")
            plt.close()
        pass


        
