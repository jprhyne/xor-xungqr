import matplotlib.pyplot as plt
import numpy as np
def lineStartsWith(line, patterns):
    if isinstance(patterns, str):
        pattern = patterns
        lenPat = len(pattern)
        if len(line) < lenPat:
            return False
        if line[0:lenPat] == pattern:
            return True
        return False
    else:
        for pattern in patterns:
            lenPat = len(pattern)
            if len(line) < lenPat:
                continue
            if line[0:lenPat] == pattern:
                return True
        return False

def createFigure(xVec, yVecs, labelVec, title, xLabel, yLabel, filename=None, xTicks=None, yTicks=None, yScale=None, nDash=0):
    fig,ax = plt.subplots()
    ax.set_title(title)
    ax.set_xlabel(xLabel)
    ax.set_ylabel(yLabel)

    numY = yVecs.shape[0]
    numL = labelVec.shape[0]

    if numY != numL:
        return

    if xTicks != None:
        ax.set_xticks(xTicks)

    if yTicks != None:
        ax.set_yticks(yTicks)
    else:
        ax.set_ylim(bottom=0, top=np.max(yVecs)*1.1)

    if yScale != None:
        ax.set_yscale(yScale)

    # Plot the first nDash elements as dashed lines
    for i in range(nDash):
        ax.plot(xVec, yVecs[i], '--', label=labelVec[i])
    # Plot the rest normally
    for i in range(nDash, numY):
        ax.plot(xVec, yVecs[i], label=labelVec[i])

    ax.legend()
    if filename != None:
        fig.savefig(filename)
