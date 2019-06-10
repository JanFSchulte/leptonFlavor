from ROOT import * 
from numpy import array as ar
from array import array
from setTDRStyle import setTDRStyle
from copy import deepcopy
import ratios
from helpersEle import *
from defsEle import getPlot, Backgrounds, Signals, Data
import uuid


def binning(channel='muon'):
    if channel == 'muon':
        nbins = 51
        m_min = 70.
        m_max = 4000.
    if channel == 'electron':
        return (range(50, 120, 5) +
                range(120, 150, 5) +
                range(150, 200, 10) +
                range(200, 600, 20) +
                range(600, 900, 30) +
                range(900, 1250, 50) +
                range(1250, 1600, 60) +
                range(1600, 1900, 70) +
                range(1900, 4000, 80) +
                range(4000, 5000, 100) +
                [5000])

    # Calculate logarithmic bins
    width = (math.log(m_max) - math.log(m_min)) / nbins
    logbins = []
    # Exceed m_max to start with Z' binning, but reach 5 TeV
    for i in range(0, nbins + 8):
        logbins.append(math.exp(math.log(m_min) + width * i))

    return logbins


def rebin(hist, binning):
    return hist.Rebin(len(binning) - 1, 'hist_' + uuid.uuid4().hex, array('d', binning))

def main():
	### for data

	
	histos = ["BBBE"]
	#~ labels = ["dielectron_Moriond2017_BB","dielectron_Moriond2017_BE"]
	labels = ["dielectron_Moriond2017_BBBE"]
	

	#~ bins = [4,7]
	bins = [4]

	massPlot = getPlot("massPlotForLimit")
	#~ massPlotSmeared = getPlot("massPlotSmeared")
	massPlotUp = getPlot("massPlotUp")
	massPlotDown = getPlot("massPlotDown")
	massPlotPUUp = getPlot("massPlotPUUp")
	massPlotPUDown = getPlot("massPlotPUDown")
	#~ massPlotWeighted = getPlot("massPlotWeighted")



	for i, histo in enumerate(histos):
		label = labels[i]
		
		drellyan = Process(getattr(Backgrounds,"DrellYan"))
		other = Process(getattr(Backgrounds,"Other"))
		#~ drellyan = Process(getattr(Signals,"CITo2Mu_Lam34TeVConLL"))

		dyHist = deepcopy(drellyan.loadHistogramProjected(massPlot, bins[i]))
		dyHistScaleUp = deepcopy(drellyan.loadHistogramProjected(massPlotUp, bins[i]))
		dyHistScaleDown = deepcopy(drellyan.loadHistogramProjected(massPlotDown, bins[i]))
		dyHistPUUp = deepcopy(drellyan.loadHistogramProjected(massPlotPUUp, bins[i]))
		dyHistPUDown = deepcopy(drellyan.loadHistogramProjected(massPlotPUDown, bins[i]))

		dyHist.Add(deepcopy(other.loadHistogramProjected(massPlot, bins[i])))
		dyHistScaleUp.Add(deepcopy(other.loadHistogramProjected(massPlotUp, bins[i])))
		dyHistScaleDown.Add(deepcopy(other.loadHistogramProjected(massPlotDown, bins[i])))
		dyHistPUUp.Add(deepcopy(other.loadHistogramProjected(massPlotPUUp, bins[i])))
		dyHistPUDown.Add(deepcopy(other.loadHistogramProjected(massPlotPUDown, bins[i])))

		dyHist.Add(deepcopy(other.loadHistogramProjected(massPlot, 7)))
		dyHistScaleUp.Add(deepcopy(other.loadHistogramProjected(massPlotUp,7)))
		dyHistScaleDown.Add(deepcopy(other.loadHistogramProjected(massPlotDown, 7)))
		dyHistPUUp.Add(deepcopy(other.loadHistogramProjected(massPlotPUUp,7)))
		dyHistPUDown.Add(deepcopy(other.loadHistogramProjected(massPlotPUDown, 7)))

		dyHist.Add(deepcopy(drellyan.loadHistogramProjected(massPlot, 7)))
		dyHistScaleUp.Add(deepcopy(drellyan.loadHistogramProjected(massPlotUp,7)))
		dyHistScaleDown.Add(deepcopy(drellyan.loadHistogramProjected(massPlotDown, 7)))
		dyHistPUUp.Add(deepcopy(drellyan.loadHistogramProjected(massPlotPUUp,7)))
		dyHistPUDown.Add(deepcopy(drellyan.loadHistogramProjected(massPlotPUDown, 7)))


		#~ rebin = 250
		
		binningEle = binning(channel="electron")
		dyHist = rebin(dyHist,binningEle)
		dyHistScaleUp = rebin(dyHistScaleUp,binningEle)
		dyHistScaleDown = rebin(dyHistScaleDown,binningEle)		
		dyHistPUUp = rebin(dyHistPUUp,binningEle)
		dyHistPUDown = rebin(dyHistPUDown,binningEle)		
		
		#~ dyHist.Rebin(rebin)
		#~ dyHistPUUp.Rebin(rebin)
		#~ dyHistPUDown.Rebin(rebin)
		#~ dyHistScaleUp.Rebin(rebin)
		#~ dyHistScaleDown.Rebin(rebin)
		hCanvas = TCanvas("hCanvas", "Distribution", 800,800)
	
		plotPad = ROOT.TPad("plotPad","plotPad",0,0.3,1,1)
		ratioPad = ROOT.TPad("ratioPad","ratioPad",0,0.,1,0.3)
		setTDRStyle()		
		plotPad.UseCurrentStyle()
		ratioPad.UseCurrentStyle()
		plotPad.Draw()	
		ratioPad.Draw()	
		plotPad.cd()
		
		plotPad.DrawFrame(250,0.00001,5000,500000000,"; dielectron mass [GeV]; Events / 250 GeV")
		plotPad.SetLogy()
		dyHist.Draw("samehist")
		dyHist.SetFillColor(kWhite)
		dyHistScaleUp.SetFillColor(kWhite)
		dyHistScaleDown.SetFillColor(kWhite)
		dyHistPUUp.SetFillColor(kWhite)
		dyHistPUDown.SetFillColor(kWhite)
		dyHistScaleUp.SetLineColor(kBlue)
		dyHistScaleDown.SetLineColor(kBlue)
		dyHistPUUp.SetLineColor(kOrange+2)
		dyHistPUDown.SetLineColor(kOrange+2)
		dyHistScaleDown.Draw("samehist")
		dyHistScaleUp.Draw("samehist")
		dyHistPUUp.Draw("samehist")
		dyHistPUDown.Draw("samehist")
		
		legend = TLegend(0.375, 0.6, 0.925, 0.925)
		legend.SetFillStyle(0)
		legend.SetBorderSize(0)
		legend.SetTextFont(42)		
		legend.AddEntry(dyHist,"Default","l")	
		legend.AddEntry(dyHistScaleDown,"Scale Uncertainty","l")	
		legend.AddEntry(dyHistPUUp,"PU Uncertainty","l")	
		legend.Draw()
		
		ratioPad.cd()
		xMin = 250
		xMax = 5000
		yMax = 1.5
		yMin = 0
		if "BE" in label:
			yMax = 1.5
			yMin = 0.5
		ratioGraphs =  ratios.RatioGraph(dyHist,dyHistScaleUp, xMin=xMin, xMax=xMax,title="Default / Uncert",yMin=yMin,yMax=yMax,ndivisions=10,color=kBlue,adaptiveBinning=10000)
		ratioGraphs.draw(ROOT.gPad,True,False,True,chi2Pos=0.8)				
		ratioGraphs2 =  ratios.RatioGraph(dyHist,dyHistScaleDown, xMin=xMin, xMax=xMax,title="Default / Uncert",yMin=0.8,yMax=1.3,ndivisions=10,color=kBlue,adaptiveBinning=10000)
		ratioGraphs2.draw(ROOT.gPad,False,False,True,chi2Pos=0.8)		
		ratioGraphs3 =  ratios.RatioGraph(dyHist,dyHistPUDown, xMin=xMin, xMax=xMax,title="Default / Uncert",yMin=0.8,yMax=1.3,ndivisions=10,color=kOrange+2,adaptiveBinning=10000)
		ratioGraphs3.draw(ROOT.gPad,False,False,True,chi2Pos=0.8)		
		ratioGraphs4 =  ratios.RatioGraph(dyHist,dyHistPUUp, xMin=xMin, xMax=xMax,title="Default / Uncert",yMin=0.8,yMax=1.3,ndivisions=10,color=kOrange+2,adaptiveBinning=10000)
		ratioGraphs4.draw(ROOT.gPad,False,False,True,chi2Pos=0.8)		
	
		#~ ratio = ratioGraphs.getGraph()

		ratio = dyHist.Clone()
		ratio.Divide(dyHistScaleUp)
		ratio.GetXaxis().SetRangeUser(0,5000)
		ratio.GetYaxis().SetRangeUser(0,2)
		func = TF1("f1","pol4",250,5000)
		#~ func.SetParameter(0,1.1)
		ratio.Fit("f1","R")
		#~ func.Draw("sameL")			
		
		
		plotF1 = TF1("plotF1","pol4",250,5000)
		plotF1.SetParameter(0,f1.GetParameter(0))
		plotF1.SetParameter(1,f1.GetParameter(1))
		plotF1.SetParameter(2,f1.GetParameter(2))
		plotF1.SetParameter(3,f1.GetParameter(3))
		plotF1.SetParameter(4,f1.GetParameter(4))
		ratio2 = dyHist.Clone()
		#~ ratio2 = ratioGraphs2.getGraph()
			
		ratio2.GetXaxis().SetRangeUser(0,5000)	
		ratio2.GetYaxis().SetRangeUser(0,2)			
		ratio2.Divide(dyHistScaleDown)
		func2 = TF1("f2","pol4",250,5000)
		#~ func2.SetParameter(0,0.9)
		ratio2.Fit("f2","R")
		#~ func2.Draw("sameL")			
		plotF2 = TF1("plotF1","pol4",250,5000)
		plotF2.SetParameter(0,f2.GetParameter(0))
		plotF2.SetParameter(1,f2.GetParameter(1))
		plotF2.SetParameter(2,f2.GetParameter(2))
		plotF2.SetParameter(3,f2.GetParameter(3))
		plotF2.SetParameter(4,f2.GetParameter(4))
		
		#~ ratioGraphs =  ratios.RatioGraph(dyHist,dyHistScaleUp, xMin=xMin, xMax=xMax,title="Default / Uncert",yMin=yMin,yMax=yMax,ndivisions=10,color=kBlue,adaptiveBinning=10000)
		ratioGraphs.draw(ROOT.gPad,True,False,True,chi2Pos=0.8)				
		#~ ratioGraphs2 =  ratios.RatioGraph(dyHist,dyHistScaleDown, xMin=xMin, xMax=xMax,title="Default / Uncert",yMin=0.8,yMax=1.3,ndivisions=10,color=kBlue,adaptiveBinning=10000)
		ratioGraphs2.draw(ROOT.gPad,False,False,True,chi2Pos=0.8)		
		#~ ratioGraphs3 =  ratios.RatioGraph(dyHist,dyHistPUDown, xMin=xMin, xMax=xMax,title="Default / Uncert",yMin=0.8,yMax=1.3,ndivisions=10,color=kOrange+2,adaptiveBinning=10000)
		ratioGraphs3.draw(ROOT.gPad,False,False,True,chi2Pos=0.8)		
		#~ ratioGraphs4 =  ratios.RatioGraph(dyHist,dyHistPUUp, xMin=xMin, xMax=xMax,title="Default / Uncert",yMin=0.8,yMax=1.3,ndivisions=10,color=kOrange+2,adaptiveBinning=10000)
		ratioGraphs4.draw(ROOT.gPad,False,False,True,chi2Pos=0.8)			
		
		plotF1.Draw("sameL")
		plotF2.Draw("sameL")
		hCanvas.Print("uncertainties_%s.pdf"%label)
		
main()
