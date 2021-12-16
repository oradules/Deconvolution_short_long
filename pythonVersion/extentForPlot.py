class extentForPlot:
    def __init__(self,f):
        delta = f[1] - f[0]
        result= [f[0] - delta/2, f[-1] + delta/2]
        self.result=result
