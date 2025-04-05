# Create hierarchy validation plots (eff, purity, completeness, vtx etc.)

import argparse
import array
import math
import os
import ROOT
import sys

class GlobalSettings:
    def __init__(self):
        # Quality cuts
        self.minCompleteness = 0.1
        self.minPurity = 0.5
        self.minNSharedHits = 5

        # List of particles
        self.pdgList = [
            (11, 'electron'),
            (13, 'muon'),
            (22, 'photon'),
            (211, 'piplus'),
            (-211, 'piminus'),
            (2212, 'proton')
        ]

        #Shared output files for overlayed histograms
        self.histMCFileName = 'Overlay_MCHistos.root'
        self.histEvtFileName = 'Overlay_EVTHistos.root'


class histoMCList(object):

    # Object storing the histograms for a given particle type
    def __init__(self, hitsAll, hitsEff, mtmAll, mtmEff, completeness, purity):
        self.hitsAll = hitsAll
        self.hitsEff = hitsEff
        self.mtmAll = mtmAll
        self.mtmEff = mtmEff
        self.completeness = completeness
        self.purity = purity
        

class histoVtxList(object):

    # Object storing the primary vertex dX, dY, dZ and dR histograms
    def __init__(self, hVtxDX, hVtxDY, hVtxDZ, hVtxDR):
        self.hVtxDX = hVtxDX
        self.hVtxDY = hVtxDY
        self.hVtxDZ = hVtxDZ
        self.hVtxDR = hVtxDR
        

class crystalBallFun(object):

    def __call__(self, xArr, pars):

        # Crystal Ball function for vertex residual fits
        x = xArr[0]
        norm = pars[0]
        x0 = pars[1]
        sigmaL = abs(pars[2])
        sigmaR = abs(pars[3])
        alphaL = abs(pars[4])
        nL = abs(pars[5])
        alphaR = abs(pars[6])
        nR = abs(pars[7])

        t = 0.0
        if x < x0:
            t = (x - x0)/sigmaL
        else:
            t = (x - x0)/sigmaR
        
        value = 0.0
        if (t < alphaL):
            value = self.getTail(t, alphaL, nL)
        elif (t <= alphaR):
            value = math.exp(-0.5*t*t)
        else:
            value = self.getTail(-t, alphaR, nR)

        return value*norm

    def getTail(self, t, alpha, n):

        result = 0.0
        if abs(alpha) > 0.0:
            a = math.pow(abs(n/alpha), n) * math.exp(-0.5*alpha*alpha)
            b = (n/alpha) - alpha
            result = a/math.pow(abs(b-t), n)
        return result



def readSampleList(filename):
    from collections import defaultdict

    sampleMap = defaultdict(dict)

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) != 3:
                print(f"Skipping malformed line: {line}")
                continue

            config, ftype, path = parts
            if ftype not in ['mc', 'evt']:
                print(f"Unknown type '{ftype}' in line: {line}")
                continue

            sampleMap[config][ftype] = path

    # Filter only complete configs (must have both mc and evt)
    sampleList = []
    for config, files in sampleMap.items():
        if 'mc' in files and 'evt' in files:
            sampleList.append((config, files['mc'], files['evt']))
        else:
            print(f"Incomplete config '{config}': missing mc or evt")

    return sampleList



def getParticleType(mcPDG):

    name = 'Unknown'
    absPDG = abs(mcPDG)
    if absPDG == 13:
        name = 'muon'
    elif absPDG == 11:
        name = 'electron'
    elif absPDG == 2212:
        name = 'proton'
    elif absPDG == 22:
        name = 'photon'
    elif mcPDG == 211:
        name = 'piplus'
    elif mcPDG == -211:
        name = 'piminus'

    return name


def defineMCHistos(settings):

    # Create empty MC histograms for each particle type.
    # Store them in the histogram map, which is returned.
    # Map key = particle type, value = histogram list object
    histMCMap = {}

    # hits binning
    nHitBins = 35
    nHitBinEdges = nHitBins + 1
    hitsBinning = array.array('d', [0.0]*nHitBinEdges)
    for iB in range(nHitBinEdges):
        edge = math.pow(10.0, 1.0 + (iB*1.0 + 2.0)*0.1)
        hitsBinning[iB] = edge
        
    # momentum binning
    nMtmBins = 26
    mtmBinning = array.array('d', [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4,
                                   1.6, 2.0, 2.4, 2.8, 3.4, 4.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0])

    # Loop over the particle types
    for iPDG,pLabel in settings.pdgList:

        # all hits, hits efficiency
        hHitsAll = ROOT.TH1F('{0}_HitsAll'.format(pLabel), '', nHitBins, hitsBinning)
        hHitsAll.SetDirectory(0)
        hHitsAll.GetXaxis().SetTitle('Number of Hits')
        hHitsAll.GetYaxis().SetTitle('Number of Events')
        
        hHitsEff = ROOT.TH1F('{0}_HitsEff'.format(pLabel), '', nHitBins, hitsBinning)
        hHitsEff.SetDirectory(0)
        hHitsEff.GetXaxis().SetTitle('Number of Hits')
        hHitsEff.GetYaxis().SetTitle('Reconstruction Efficiency')

        # all momentum, momentum efficiency
        hMtmAll = ROOT.TH1F('{0}_MtmAll'.format(pLabel), '', nMtmBins, mtmBinning)
        hMtmAll.SetDirectory(0)
        hMtmAll.GetXaxis().SetTitle('True Momentum [GeV]')
        hMtmAll.GetYaxis().SetTitle('Number of Events')

        hMtmEff = ROOT.TH1F('{0}_MtmEff'.format(pLabel), '', nMtmBins, mtmBinning)
        hMtmEff.SetDirectory(0)
        hMtmEff.GetXaxis().SetTitle('True Momentum [GeV]')
        hMtmEff.GetYaxis().SetTitle('Reconstruction Efficiency')

        # Completeness
        hCompleteness = ROOT.TH1F('{0}_Completeness'.format(pLabel), '', 51, -0.01, 1.01)
        hCompleteness.SetDirectory(0)
        hCompleteness.GetXaxis().SetTitle('Completeness')
        hCompleteness.GetYaxis().SetTitle('Fraction of Events')

        # Purity
        hPurity = ROOT.TH1F('{0}_Purity'.format(pLabel), '', 51, -0.01, 1.01)
        hPurity.SetDirectory(0)
        hPurity.GetXaxis().SetTitle('Purity')
        hPurity.GetYaxis().SetTitle('Fraction of Events')

        histMCMap[pLabel] = histoMCList(hHitsAll, hHitsEff, hMtmAll, hMtmEff,
                                        hCompleteness, hPurity)

    # Return the histogram map
    return histMCMap


def defineVtxHistos():

    # Define primary vertex histograms for all interactions
    hVtxDX = ROOT.TH1F('allVtxDX', '', 100, -5.0, 5.0)
    hVtxDX.SetDirectory(0)
    hVtxDX.GetXaxis().SetTitle('Vertex #DeltaX [cm]')
    hVtxDX.GetYaxis().SetTitle('Number of Events')

    hVtxDY = ROOT.TH1F('allVtxDY', '', 100, -5.0, 5.0)
    hVtxDY.SetDirectory(0)
    hVtxDY.GetXaxis().SetTitle('Vertex #DeltaY [cm]')
    hVtxDY.GetYaxis().SetTitle('Number of Events')

    hVtxDZ = ROOT.TH1F('allVtxDZ', '', 100, -5.0, 5.0)
    hVtxDZ.SetDirectory(0)
    hVtxDZ.GetXaxis().SetTitle('Vertex #DeltaZ [cm]')
    hVtxDZ.GetYaxis().SetTitle('Number of Events')

    hVtxDR = ROOT.TH1F('allVtxDR', '', 100, 0.0, 5.0)
    hVtxDR.SetDirectory(0)
    hVtxDR.GetXaxis().SetTitle('Vertex #DeltaR [cm]')
    hVtxDR.GetYaxis().SetTitle('Number of Events')

    # Object storing the list of histograms
    hVtxList = histoVtxList(hVtxDX, hVtxDY, hVtxDZ, hVtxDR)
    return hVtxList
    

def createMCHistos(pars):

    print('mcFileName = {0}, mcTreeName = {1}'.format(pars.mcFileName, pars.mcTreeName))

    # Define performance histograms for each particle type
    histMCMap = defineMCHistos(pars)

    # Open MC hierachy file and its tree
    mcFile = ROOT.TFile.Open(pars.mcFileName, 'read')
    mcTree = mcFile.Get(pars.mcTreeName)

    # Loop over MC hierarchy tree entries
    nMC = mcTree.GetEntries()
    nTotPass = 0
    nTotFail = 0
    nTotal = 0

    for i in range(nMC):

        # Get tree entry
        mcTree.GetEntry(i)

        event = getattr(mcTree, 'event')

        if i%10000 == 0:
            print('MC entry {0}'.format(nMC-i))
            
        # Get MC particle PDG Id
        mcPDG = getattr(mcTree, 'mcPDG')
        # Get particle name type
        mcType = getParticleType(mcPDG)

        # Require known particle type
        if mcType == 'Unknown':
            continue

        # Histogram list for given particle type
        hList = histMCMap[mcType]

        # Number of hits and true momentum
        hHitsAll = hList.hitsAll
        mcNHits = getattr(mcTree, 'mcNHits')
        hHitsAll.Fill(mcNHits)

        mcMtm = getattr(mcTree, 'mcMtm')
        nMatches = getattr(mcTree, 'nMatches')
        
        hMtmAll = hList.mtmAll
        hMtmAll.Fill(mcMtm)

        # Find best completeness MC-reco match and its corresponding purity
        completeVect = getattr(mcTree, 'completenessVector')
        purityVect = getattr(mcTree, 'purityVector')
        nSharedHitsVect = getattr(mcTree, 'nSharedHitsVector')
            
        complete = 0.0
        purity = 0.0
        nSharedHits = 0
        bestNShared = 0
        for iC,nShared in enumerate(nSharedHitsVect):
            if nShared > bestNShared:
                bestNShared = nShared
                complete = completeVect[iC]
                purity = purityVect[iC]

        # Check for minimum completeness, purity and number of shared hits
        if nMatches > 0 and complete >= pars.minCompleteness and purity >= pars.minPurity \
           and bestNShared >= pars.minNSharedHits:

            # Efficiency numerator (hits & true momentum)
            hHitsEff = hList.hitsEff
            hHitsEff.Fill(mcNHits)

            hMtmEff = hList.mtmEff
            hMtmEff.Fill(mcMtm)

            # Fill completeness and purity histos
            hCompleteness = hList.completeness
            hCompleteness.Fill(complete)

            hPurity = hList.purity
            hPurity.Fill(purity)


    # Process hits and momentum histograms to get their efficiencies,
    # and normalise the completeness and purity histograms.
    # Then, write these all to the histogram output file
    print('Creating {0}'.format(pars.histMCFileName))
    hMCOutFile = ROOT.TFile.Open(pars.histMCFileName, 'recreate')

    # Loop over particle types
    for pdgId, pLabel in pars.pdgList:

        # Histogram list for given particle type
        hList = histMCMap[pLabel]

        # Find efficiencies
        setEffHist(hList.hitsEff, hList.hitsAll)
        setEffHist(hList.mtmEff, hList.mtmAll)

        # Normalise completeness and purity
        nCompleteness = hList.completeness.GetEntries()*1.0
        if nCompleteness > 0.0:
            hList.completeness.Scale(1.0/nCompleteness)

        nPurity = hList.purity.GetEntries()*1.0
        if nPurity > 0.0:
            hList.purity.Scale(1.0/nPurity)

        # Write out the histograms
        hMCOutFile.cd()
        hList.hitsAll.Write()
        hList.hitsEff.Write()
        hList.mtmAll.Write()
        hList.mtmEff.Write()
        hList.completeness.Write()
        hList.purity.Write()
        
    # Close the files
    hMCOutFile.Close()
    mcFile.Close()


def createVtxHistos(pars):
    
    # Event vertex histograms
    hVtxList = defineVtxHistos(pars)

    # Open event hierarchy file
    evtFile = ROOT.TFile.Open(pars.evtFileName, 'read')
    evtTree = evtFile.Get(pars.evtTreeName)

    # Loop over event entries
    nEvt = evtTree.GetEntries()
    for i in range(nEvt):

        # Get tree entry
        evtTree.GetEntry(i)

        if i%10000 == 0:
            print('Evt entry {0}'.format(nEvt-i))

        # Get vertex residuals
        vtxDx = getattr(evtTree, 'vtxDx')
        vtxDy = getattr(evtTree, 'vtxDy')
        vtxDz = getattr(evtTree, 'vtxDz')
        vtxDr = getattr(evtTree, 'vtxDr')

        hVtxList.hVtxDX.Fill(vtxDx)
        hVtxList.hVtxDY.Fill(vtxDy)
        hVtxList.hVtxDZ.Fill(vtxDz)
        hVtxList.hVtxDR.Fill(vtxDr)

    # Write out histograms
    print('Creating {0}'.format(pars.histEvtFileName))
    hEvtOutFile = ROOT.TFile.Open(pars.histEvtFileName, 'recreate')
    hEvtOutFile.cd()
    hVtxList.hVtxDX.Write()
    hVtxList.hVtxDY.Write()
    hVtxList.hVtxDZ.Write()
    hVtxList.hVtxDR.Write()

    # Close files
    hEvtOutFile.Close()
    evtFile.Close()
        
            
def setEffHist(hEff, hAll):

    # Modify the numerator hEff histogram to be the efficiency
    # by dividing by the denominator hAll histogram bin contents.
    # Assumes both have the same binning, which they should do

    # Loop over bins, including under and overflow
    nBins = hEff.GetXaxis().GetNbins()
    for i in range(-1, nBins+1):
        i1 = i + 1
        num = hEff.GetBinContent(i1)
        denom = hAll.GetBinContent(i1)
        # Efficiency and its binomial error
        eff = (num/denom) if (denom > 0.0) else 0.0
        err = math.sqrt(eff*(1.0 - eff)/denom) if (denom > 0.0) else 0.0
        # Update the efficiency bin content
        hEff.SetBinContent(i1, eff)
        hEff.SetBinError(i1, err)


def setCBFun(hist, cbFun, theMean = -999.0, theSigma = -999.0):

    # Set the Crystal Ball fit function for the given histogram.
    # Note that we pass cbFun = crystalBallFun() as a callable argument
    # so that it remains in scope when the function is returned
    histName = hist.GetName()
    funName = '{0}Fun'.format(histName)
    xAxis = hist.GetXaxis()
    xMin = xAxis.GetXmin()
    xMax = xAxis.GetXmax()
    vtxDXFun = ROOT.TF1(funName, cbFun, xMin, xMax, 8)

    mean = theMean if (theMean > -999.0) else hist.GetMean()
    sigma = theSigma if (theSigma > -999.0) else hist.GetStdDev()
    norm = hist.GetMaximum()

    fun = ROOT.TF1(funName, cbFun, xMin, xMax, 8)
    fun.SetParameters(norm, mean, sigma, sigma, 1.2, 1.0, 1.2, 1.0)
    fun.SetParNames('N', '#mu', '#sigma_{L}', '#sigma_{R}',
                    '#alpha_{L}', 'n_{L}', '#alpha_{R}', 'n_{R}')
    
    return fun

def fillVtxHistosFromTree(evtFileName, evtTreeName, vtxHists):
    print(f"→ Reading event file: {evtFileName}")
    evtFile = ROOT.TFile.Open(evtFileName, "read")
    evtTree = evtFile.Get(evtTreeName)

    if not evtTree:
        print(f"Could not find tree '{evtTreeName}' in file {evtFileName}")
        return

    nEvt = evtTree.GetEntries()
    print(f"  Entries: {nEvt}")

    for i in range(nEvt):
        evtTree.GetEntry(i)

        if i % 10000 == 0:
            print(f"  Remaining entries: {nEvt - i}")

        dx = getattr(evtTree, 'vtxDx')
        dy = getattr(evtTree, 'vtxDy')
        dz = getattr(evtTree, 'vtxDz')
        dr = getattr(evtTree, 'vtxDr')

        vtxHists.hVtxDX.Fill(dx)
        vtxHists.hVtxDY.Fill(dy)
        vtxHists.hVtxDZ.Fill(dz)
        vtxHists.hVtxDR.Fill(dr)

    evtFile.Close()

def fillMCHistosFromTree(mcFileName, mcTreeName, histMCMap, settings):
    print(f'→ Reading MC file: {mcFileName}')
    mcFile = ROOT.TFile.Open(mcFileName, 'read')
    mcTree = mcFile.Get(mcTreeName)

    if not mcTree:
        print(f"Could not find tree '{mcTreeName}' in file {mcFileName}")
        return

    nMC = mcTree.GetEntries()
    print(f'  Entries: {nMC}')

    for i in range(nMC):
        mcTree.GetEntry(i)

        if i % 10000 == 0:
            print(f'  Remaining entries: {nMC - i}')

        mcPDG = getattr(mcTree, 'mcPDG')
        mcType = getParticleType(mcPDG)
        if mcType == 'Unknown' or mcType not in histMCMap:
            continue

        hList = histMCMap[mcType]

        mcNHits = getattr(mcTree, 'mcNHits')
        mcMtm = getattr(mcTree, 'mcMtm')
        nMatches = getattr(mcTree, 'nMatches')

        hList.hitsAll.Fill(mcNHits)
        hList.mtmAll.Fill(mcMtm)

        completeVect = getattr(mcTree, 'completenessVector')
        purityVect = getattr(mcTree, 'purityVector')
        nSharedHitsVect = getattr(mcTree, 'nSharedHitsVector')

        # Best match
        bestNShared = 0
        bestComplete = 0.0
        bestPurity = 0.0

        for iC, nShared in enumerate(nSharedHitsVect):
            if nShared > bestNShared:
                bestNShared = nShared
                bestComplete = completeVect[iC]
                bestPurity = purityVect[iC]

        # Apply quality cuts
        if (nMatches > 0 and
            bestComplete >= settings.minCompleteness and
            bestPurity >= settings.minPurity and
            bestNShared >= settings.minNSharedHits):

            hList.hitsEff.Fill(mcNHits)
            hList.mtmEff.Fill(mcMtm)
            hList.completeness.Fill(bestComplete)
            hList.purity.Fill(bestPurity)

    mcFile.Close()

    # Finalize histograms: normalize and compute efficiency
    for pdgId, pLabel in settings.pdgList:
        if pLabel not in histMCMap:
            continue

        hList = histMCMap[pLabel]

        # Compute efficiencies
        setEffHist(hList.hitsEff, hList.hitsAll)
        setEffHist(hList.mtmEff, hList.mtmAll)

        # Normalize completeness and purity
        nComp = hList.completeness.GetEntries()
        if nComp > 0:
            hList.completeness.Scale(1.0 / nComp)

        nPur = hList.purity.GetEntries()
        if nPur > 0:
            hList.purity.Scale(1.0 / nPur)

def writeVtxHistos(vtxHists, outputFileName):
    print(f"→ Writing vertex histograms to: {outputFileName}")
    outFile = ROOT.TFile.Open(outputFileName, "recreate")

    outFile.cd()
    vtxHists.hVtxDX.Write()
    vtxHists.hVtxDY.Write()
    vtxHists.hVtxDZ.Write()
    vtxHists.hVtxDR.Write()

    outFile.Close()

def writeMcHistos(histMCMap, outputFileName, settings):
    print(f"→ Writing histograms to: {outputFileName}")
    outFile = ROOT.TFile.Open(outputFileName, 'recreate')

    for pdgId, pLabel in settings.pdgList:
        if pLabel not in histMCMap:
            continue
        hList = histMCMap[pLabel]

        outFile.cd()
        hList.hitsAll.Write()
        hList.hitsEff.Write()
        hList.mtmAll.Write()
        hList.mtmEff.Write()
        hList.completeness.Write()
        hList.purity.Write()

    outFile.Close()
    print(f"Done writing to: {outputFileName}")


def run(args):

    settings = GlobalSettings()
    sampleList = readSampleList(args.sampleList)
    print("Parsed sample list:")
    for label, mcFile, evtFile in sampleList:
        print(f"  {label}:\n    MC : {mcFile}\n    EVT: {evtFile}")

    # Create the histograms. Write them to the ROOT output file with the name
    # "mcFileName_Histos.root", where mcFileName has the .root extension removed
    if args.createHistos == 1:

        allMCConfigs = {}

        # Step 1: define shared histograms
        mcHists = defineMCHistos(settings)
        vtxHists = defineVtxHistos()

        configNames = []

        for config, mcFile, evtFile in sampleList:

            print(f"→ Processing config: {config}")
            fillMCHistosFromTree(mcFile, args.mcTreeName, mcHists, settings)
            writeMcHistos(mcHists, f"{config}_MCHistos.root", settings)

            fillVtxHistosFromTree(evtFile, args.evtTreeName, vtxHists)
            writeVtxHistos(vtxHists, f"{config}_VtxHistos.root")
            
            configNames.append(config)

    
def processArgs(parser):

    # Process script arguments
    parser.add_argument('--mcTreeName', default='MC', metavar='treeName',
                        help='MC hierarchy ROOT tree [default "MC"]')

    parser.add_argument('--evtTreeName', default='Events', metavar='treeName',
                        help='Event hierarchy ROOT tree [default "Events"]')

    parser.add_argument('--createHistos', default=1, metavar='int', type=int,
                        help='Recreate histograms [1 = Yes (default), 0 = No]')

    parser.add_argument('--sampleList', default='samples.txt', metavar='file',
                        help='Text file listing sample names and .root paths')

if __name__ == '__main__':

    # Process the command line arguments
    # Use "python hierarchyPlots.py --help" to see the full list
    parser = argparse.ArgumentParser(description='List of arguments')
    processArgs(parser)
    args = parser.parse_args()

    run(args)
