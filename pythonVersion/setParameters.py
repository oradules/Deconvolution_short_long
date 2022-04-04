import numpy as np

class setParameters:
    def __init__(self, parameterFilePath,parameterFileName):
        param = np.zeros(10)
        param[6]=1
        for i in range(6):
            if i==0:
                print('Enter the Polymerase speed')
                param[i] = input() #Polym_speed
            if i==1:
                print('Enter the length of the First part of polymerase signal')
                param[i] = input() #TaillePreMarq
            if i==2:
                print('Enter the length of the Second part of polymerase signal containing the sequence for MS2 loops')
                param[i] = input() #TailleSeqMarq
            if i==3:
                print('Enter the Last part of polymerase signal where the PolII waits')
                param[i] = input() #TaillePostMarq
            if i==4:
                print('Enter the minimum inter POLII space in bp')
                param[i] = input() #EspaceInterPolyMin
            if i==5:
                print('Enter the frame length in seconds')
                param[i] = input()
        param[7] = 1/param[5] # FreqEchImg = (1/FrameLen) 
        param[8] = (param[1] + param[2] + param[3])/param[0]# DureeSignal = (TaillePreMarq + TailleSeqMarq + TaillePostMarq) / Polym_speed
        param[9] = 1/(param[4]/param[0]) # FreqEchSimu = 1/(EspaceInterPolyMin/Polym_speed)
        np.savez(parameterFilePath+parameterFileName, 
                 Polym_speed = param[0], 
                 TaillePreMarq = param[1],
                 TailleSeqMarq = param[2],
                 TaillePostMarq = param[3],
                 EspaceInterPolyMin = param[4],
                 FrameLen = param[5],
                 Intensity_for_1_Polym = param[6],
                 FreqEchImg = param[7],
                 DureeSignal = param[8],
                 FreqEchSimu = param[9]
                )
        self.Polym_speed = param[0]
        self.TaillePreMarq = param[1]
        self.TailleSeqMarq = param[2]
        self.TaillePostMarq = param[3]
        self.EspaceInterPolyMin = param[4]
        self.FrameLen = param[5]
        self.Intensity_for_1_Polym = param[6]
        self.FreqEchImg = param[7]
        self.DureeSignal = param[8]
        self.FreqEchSimu = param[9]
