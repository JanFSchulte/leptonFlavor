from ROOT import * 
from numpy import array as ar
from array import array
from setTDRStyle import setTDRStyle
from copy import deepcopy
import ratios
from helpers import *
from defs import getPlot, Backgrounds, Signals, Data

def main():
	### for data

	
	histos = ["BB","BE"]
	labels = ["dimuon_Moriond2017_BB","dimuon_Moriond2017_BE"]
	

	bins = [4,7]

	massPlot = getPlot("massPlotForLimit")
	massPlotSmeared = getPlot("massPlotSmeared")
	massPlotUp = getPlot("massPlotUp")
	massPlotDown = getPlot("massPlotDown")
	massPlotPUUp = getPlot("massPlotPUUp")
	massPlotPUDown = getPlot("massPlotPUDown")
	massPlotWeighted = getPlot("massPlotWeighted")



	for i, histo in enumerate(histos):
		label = labels[i]
		
		drellyan = Process(getattr(Backgrounds,"DrellYan"))
		other = Process(getattr(Backgrounds,"Other"))
		#~ drellyan = Process(getattr(Signals,"CITo2Mu_Lam34TeVConLL"))

		dyHist = deepcopy(drellyan.loadHistogramProjected(massPlot, bins[i]))
		dyHistSmear = deepcopy(drellyan.loadHistogramProjected(massPlotSmeared, bins[i]))
		dyHistScaleUp = deepcopy(drellyan.loadHistogramProjected(massPlotUp, bins[i]))
		dyHistScaleDown = deepcopy(drellyan.loadHistogramProjected(massPlotDown, bins[i]))
		dyHistPUUp = deepcopy(drellyan.loadHistogramProjected(massPlotPUUp, bins[i]))
		dyHistPUDown = deepcopy(drellyan.loadHistogramProjected(massPlotPUDown, bins[i]))
		dyHistWeighted = deepcopy(drellyan.loadHistogramProjected(massPlotWeighted, bins[i]))
		if "_BE" in label:
			dyHist.Add(deepcopy(drellyan.loadHistogramProjected(massPlot, 10)))
			dyHistSmear.Add(deepcopy(drellyan.loadHistogramProjected(massPlotSmeared, 10)))
			dyHistScaleUp.Add(deepcopy(drellyan.loadHistogramProjected(massPlotUp, 10)))
			dyHistScaleDown.Add(deepcopy(drellyan.loadHistogramProjected(massPlotDown, 10)))
			dyHistPUUp.Add(deepcopy(drellyan.loadHistogramProjected(massPlotPUUp, 10)))
			dyHistPUDown.Add(deepcopy(drellyan.loadHistogramProjected(massPlotPUDown, 10)))
			dyHistWeighted.Add(deepcopy(drellyan.loadHistogramProjected(massPlotWeighted, 10)))
		dyHist.Add(deepcopy(other.loadHistogramProjected(massPlot, bins[i])))
		dyHistSmear.Add(deepcopy(other.loadHistogramProjected(massPlotSmeared, bins[i])))
		dyHistScaleUp.Add(deepcopy(other.loadHistogramProjected(massPlotUp, bins[i])))
		dyHistScaleDown.Add(deepcopy(other.loadHistogramProjected(massPlotDown, bins[i])))
		dyHistPUUp.Add(deepcopy(other.loadHistogramProjected(massPlotPUUp, bins[i])))
		dyHistPUDown.Add(deepcopy(other.loadHistogramProjected(massPlotPUDown, bins[i])))
		dyHistWeighted.Add(deepcopy(other.loadHistogramProjected(massPlotWeighted, bins[i])))
		if "_BE" in label:
			dyHist.Add(deepcopy(other.loadHistogramProjected(massPlot, 10)))
			dyHistSmear.Add(deepcopy(other.loadHistogramProjected(massPlotSmeared, 10)))
			dyHistScaleUp.Add(deepcopy(other.loadHistogramProjected(massPlotUp, 10)))
			dyHistScaleDown.Add(deepcopy(other.loadHistogramProjected(massPlotDown, 10)))
			dyHistPUUp.Add(deepcopy(other.loadHistogramProjected(massPlotPUUp, 10)))
			dyHistPUDown.Add(deepcopy(other.loadHistogramProjected(massPlotPUDown, 10)))
			dyHistWeighted.Add(deepcopy(other.loadHistogramProjected(massPlotWeighted, 10)))
		rebin = 500
		dyHist.Rebin(rebin)
		dyHistSmear.Rebin(rebin)
		dyHistPUUp.Rebin(rebin)
		dyHistPUDown.Rebin(rebin)
		dyHistScaleDown.Rebin(rebin)
		dyHistWeighted.Rebin(rebin)
		hCanvas = TCanvas("hCanvas", "Distribution", 800,800)
	
		plotPad = ROOT.TPad("plotPad","plotPad",0,0.3,1,1)
		ratioPad = ROOT.TPad("ratioPad","ratioPad",0,0.,1,0.3)
		setTDRStyle()		
		plotPad.UseCurrentStyle()
		ratioPad.UseCurrentStyle()
		plotPad.Draw()	
		ratioPad.Draw()	
		plotPad.cd()
		
		plotPad.DrawFrame(0,0.0001,5000,500000,"; dimuon mass [GeV]; Events / 500 GeV")
		plotPad.SetLogy()
		dyHist.Draw("samehist")
		dyHist.SetFillColor(kWhite)
		dyHistSmear.SetFillColor(kWhite)
		dyHistScaleDown.SetFillColor(kWhite)
		dyHistWeighted.SetFillColor(kWhite)
		dyHistPUUp.SetFillColor(kWhite)
		dyHistPUDown.SetFillColor(kWhite)
		dyHistSmear.SetLineColor(kRed)
		dyHistScaleDown.SetLineColor(kBlue)
		dyHistWeighted.SetLineColor(kGreen+2)
		dyHistPUUp.SetLineColor(kOrange+2)
		dyHistSmear.Draw("samehist")
		dyHistScaleDown.Draw("samehist")
		dyHistWeighted.Draw("samehist")
		dyHistPUUp.Draw("samehist")
		
		legend = TLegend(0.375, 0.6, 0.925, 0.925)
		legend.SetFillStyle(0)
		legend.SetBorderSize(0)
		legend.SetTextFont(42)		
		legend.AddEntry(dyHist,"Default","l")	
		legend.AddEntry(dyHistSmear,"Resolution Uncertainty","l")	
		legend.AddEntry(dyHistScaleDown,"Scale Uncertainty","l")	
		legend.AddEntry(dyHistWeighted,"ID Uncertainty","l")	
		legend.AddEntry(dyHistPUUp,"PU Uncertainty","l")	
		legend.Draw()
		
		ratioPad.cd()
		xMin = 0
		xMax = 5000
		yMax = 1.1
		yMin = 0.9
		if "BE" in label:
			yMax = 1.2
			yMin = 0.8
		ratioGraphs =  ratios.RatioGraph(dyHist,dyHistSmear, xMin=xMin, xMax=xMax,title="Default / Uncert",yMin=yMin,yMax=yMax,ndivisions=10,color=kRed,adaptiveBinning=10000)
		ratioGraphs.draw(ROOT.gPad,True,False,True,chi2Pos=0.8)		
		ratioGraphs2 =  ratios.RatioGraph(dyHist,dyHistScaleDown, xMin=xMin, xMax=xMax,title="Default / Uncert",yMin=0.8,yMax=1.3,ndivisions=10,color=kBlue,adaptiveBinning=10000)
		ratioGraphs2.draw(ROOT.gPad,False,False,True,chi2Pos=0.8)
		ratioGraphs3 =  ratios.RatioGraph(dyHist,dyHistWeighted, xMin=xMin, xMax=xMax,title="Default / Uncert",yMin=0.8,yMax=1.3,ndivisions=10,color=kGreen+2,adaptiveBinning=10000)
		ratioGraphs3.draw(ROOT.gPad,False,False,True,chi2Pos=0.8)		
		ratioGraphs4 =  ratios.RatioGraph(dyHist,dyHistPUUp, xMin=xMin, xMax=xMax,title="Default / Uncert",yMin=0.8,yMax=1.3,ndivisions=10,color=kOrange+2,adaptiveBinning=10000)
		ratioGraphs4.draw(ROOT.gPad,False,False,True,chi2Pos=0.8)		
		
		hCanvas.Print("uncertainties_%s.pdf"%label)
		
main()
