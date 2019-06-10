from ROOT import * 
from numpy import array as ar
from array import array
from setTDRStyle import setTDRStyle
from copy import deepcopy
import ratios
from helpersEle import *
from defsEle import getPlot, Backgrounds, Signals, Data

def main():
	### for data

	
	histos = ["BB","BE"]
	labels = ["dielectron_Moriond2017_BB","dielectron_Moriond2017_BE"]
	#~ labels = ["dielectron_Moriond2017_BE"]
	

	bins = [4,7]
	#~ bins = [7]

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


		rebin = 250
		dyHist.Rebin(rebin)
		dyHistPUUp.Rebin(rebin)
		dyHistPUDown.Rebin(rebin)
		dyHistScaleUp.Rebin(rebin)
		dyHistScaleDown.Rebin(rebin)
		hCanvas = TCanvas("hCanvas", "Distribution", 800,800)
	
		plotPad = ROOT.TPad("plotPad","plotPad",0,0.3,1,1)
		ratioPad = ROOT.TPad("ratioPad","ratioPad",0,0.,1,0.3)
		setTDRStyle()		
		plotPad.UseCurrentStyle()
		ratioPad.UseCurrentStyle()
		plotPad.Draw()	
		ratioPad.Draw()	
		plotPad.cd()
		
		plotPad.DrawFrame(0,0.00001,5000,500000000,"; dielectron mass [GeV]; Events / 250 GeV")
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
		xMin = 0
		xMax = 5000
		yMax = 1.5
		yMin = 0.5
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
	

				
		
		#~ plotF1.Draw("sameL")
		#~ plotF2.Draw("sameL")
		hCanvas.Print("uncertainties_%s.pdf"%label)
		
main()
