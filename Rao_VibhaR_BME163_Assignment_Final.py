import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
import argparse

plt.style.use('BME163')

def readPsl(file):
    dataList = []
    with open(file) as f:
        for line in f:
            splitLine = line.strip().split('\t')
            readChromosome, readStart, readEnd = splitLine[13], int(splitLine[15]), int(splitLine[16])
            if chromosome == readChromosome:
                keep = False
                if start < readStart < end:
                    keep = True
                if start < readEnd < end:
                    keep = True
                if readStart < start and readEnd > end:
                    keep = True
                if keep:
                    blockStarts = np.array(splitLine[20].split(',')[:-1], dtype=int)
                    blockWidths = np.array(splitLine[18].split(',')[:-1], dtype=int)
                    dataList.append([readStart, readEnd, list(blockStarts), list(blockWidths), ['exon'] * len(list(blockWidths)), 0])
    return dataList

def readGtf(file):
    gtfDict = {}
    with open(file) as f:
        for line in f:
            if not line.startswith('#'):
                splitLine = line.strip().split('\t')
                if splitLine[2] == 'exon' or splitLine[2] == 'CDS':
                    Gchromosome = splitLine[0]
                    Gstart, Gend, feature = int(splitLine[3]), int(splitLine[4]), splitLine[2]
                    info = splitLine[8]
                    transcript_id = info.split('transcript_id ')[1].split(';')[0].strip('"')
                    if transcript_id not in gtfDict:
                        gtfDict[transcript_id] = []
                    gtfDict[transcript_id].append((Gchromosome, Gstart, Gend, feature))

    gtfList = []
    for transcript_id, features in gtfDict.items():
        startsNends = []
        blockStarts = []
        blockWidths = []
        blockTypes = []
        for Gchromosome, Gstart, Gend, type1 in features:
            startsNends.append(Gstart)
            startsNends.append(Gend)
            blockStarts.append(Gstart)
            blockWidths.append(Gend - Gstart)
            blockTypes.append(type1)

        readStart = min(startsNends)
        readEnd = max(startsNends)
        if chromosome == Gchromosome:
            keep = False
            if start < readStart < end:
                keep = True
            if start < readEnd < end:
                keep = True
            if keep:
                gtfList.append([readStart, readEnd, blockStarts, blockWidths, blockTypes, 0])
    return gtfList

def stackStuff(dataList):
    rows = []

    for read in dataList:
        assigned = False
        for rowIndex in range(len(rows)):
            row = rows[rowIndex]
            previous_end = -1
            for r in row:
                previous_end = max(previous_end, r[1])
            if read[0] > previous_end:
                row.append(read)
                read[5] = rowIndex
                assigned = True
                break
        if not assigned:
            rows.append([read])
            read[5] = len(rows) - 1
    return dataList

def plotStuffGtf(dataList, panel, color, linewidth):
    dataList = stackStuff(dataList)
    for readStart, readEnd, blockStarts, blockWidths, blockTypes, yPos in dataList:
        maxHeight = max([0.25 if blockType == 'exon' else 0.5 for blockType in blockTypes])
        rectangle = mplpatches.Rectangle([readStart, yPos + 0.2], readEnd - readStart, 0.05, facecolor=color, edgecolor='black', linewidth=linewidth)
        panel.add_patch(rectangle)
        for i in range(len(blockStarts)):
            blockStart = blockStarts[i]
            blockWidth = blockWidths[i]
            height = 0.25 if blockTypes[i] == 'exon' else 0.5
            centerOffset = (0.5 - height) / 2
            rectangle = mplpatches.Rectangle([blockStart, yPos+centerOffset], blockWidth, height, facecolor=color, edgecolor='black', linewidth=linewidth)
            panel.add_patch(rectangle)

def plotStuff(dataList, panel, color, linewidth, sortBy):
    dataList = sorted(dataList, key=lambda x: x[1] if sortBy == 'end' else x[0])
    dataList = stackStuff(dataList)

    for readStart, readEnd, blockStarts, blockWidths, blockTypes, yPos in dataList:
        rectangle = mplpatches.Rectangle([readStart, yPos + 0.2], readEnd - readStart, 0.05, facecolor=color, linewidth=0)
        panel.add_patch(rectangle)
        for i in range(len(blockStarts)):
            blockStart = blockStarts[i]
            blockWidth = blockWidths[i]
            height = 0.5  
            rectangle = mplpatches.Rectangle([blockStart, yPos], blockWidth, height, facecolor=color, linewidth=0)
            panel.add_patch(rectangle)

def calculateCoverage(dataList, start, end):
    coverage = np.zeros(end - start + 1)
    for readStart, readEnd, blockStarts, blockWidths, blockTypes, yPos in dataList:
        for i in range(len(blockStarts)):
            if blockTypes[i] == 'exon': 
                blockStart = blockStarts[i]
                blockEnd = blockStart + blockWidths[i]
                for pos in range(max(blockStart, start), min(blockEnd, end)):
                    coverage[pos - start] += 1
    return coverage

def plotCoverageHistogram(panel, start, end, coverage, color):
    bins = np.arange(start, end + 1)
    hist, _ = np.histogram(bins, bins=bins, weights=coverage)
    max_height = max(hist)
    for i in range(len(hist)):
        height = hist[i]
        if height > 0:  
            rect = mplpatches.Rectangle((start + i, panel.get_ylim()[1] - height), 1, height, facecolor=color, edgecolor=color, linewidth=0)
            panel.add_patch(rect)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p5',type=str, action='store')
    parser.add_argument('-p6',type=str, action='store')
    parser.add_argument('-g',type=str, action='store')
    parser.add_argument('-c',type=str, action='store')
    parser.add_argument('-o',type=str, action='store',default='Rao_VibhaR_BME163_Assignment_Final.png')
    args = parser.parse_args()

    global chromosome, start, end
    coord = args.c.split(':')
    chromosome = coord[0]
    start, end = map(int, coord[1].split('-'))

    figureWidth = 5
    figureHeight = 6
    plt.figure(figsize=(figureWidth, figureHeight))

    panel0 = plt.axes([0.1 / figureWidth, 0.1 / figureHeight, 4 / figureWidth, 0.4 / figureHeight], frameon=True)
    panel1 = plt.axes([0.1 / figureWidth, 0.5 / figureHeight, 4 / figureWidth, 1.5 / figureHeight], frameon=True)
    panel2 = plt.axes([0.1 / figureWidth, 2.2 / figureHeight, 4 / figureWidth, 1.5 / figureHeight], frameon=True)
    panel3 = plt.axes([0.1 / figureWidth, 3.9 / figureHeight, 4 / figureWidth, 1.5 / figureHeight], frameon=True)

    panel0.tick_params(bottom=False, labelbottom=False, left=False, labelleft=False,
                       right=False, labelright=False, top=False, labeltop=False)
    panel1.tick_params(bottom=False, labelbottom=False, left=False, labelleft=False,
                       right=False, labelright=False, top=False, labeltop=False)
    panel2.tick_params(bottom=False, labelbottom=False, left=False, labelleft=False,
                       right=False, labelright=False, top=False, labeltop=False)
    panel3.tick_params(bottom=False, labelbottom=False, left=False, labelleft=False,
                       right=False, labelright=False, top=False, labeltop=False)

    psl5 = readPsl(args.p5)
    psl6 = readPsl(args.p6)
    gtf = readGtf(args.g)

    iblue = (88 / 255, 85 / 255, 120 / 255)
    iorange = (230 / 255, 87 / 255, 43 / 255)
    
    plotStuffGtf(gtf, panel3, 'grey', 0.25)
    plotStuff(psl5, panel2, iorange, 0.05, 'end')
    plotStuff(psl6, panel1, iblue, 0.05, 'start')

    coverage = calculateCoverage(psl6, start, end)
    panel0.set_ylim(0, max(coverage)) 
    plotCoverageHistogram(panel0, start, end, coverage, iblue)

    panel3.set_xlim(start - (end - start) * 0.01, end + (end - start) * 0.01)  # Adding padding for x-axis
    panel2.set_xlim(start, end)
    panel1.set_xlim(start, end)
    panel0.set_xlim(start, end)

    panel3.set_ylim(-0.75, (max([read[5] for read in gtf]) + 1.25) + max([read[5] for read in gtf]) * 0.1)
    panel2.set_ylim(0.25, (max([read[5] for read in psl5]) + 1.25) + max([read[5] for read in psl5]) * 0.1)
    panel1.set_ylim(0.25, (max([read[5] for read in psl6]) + 1.25) + max([read[5] for read in psl6]) * 0.1)
     

    plt.savefig(args.o, dpi=2400)

if __name__ == "__main__":
    main()
