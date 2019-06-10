from ROOT import * 
from numpy import array as ar
from array import array
from setTDRStyle import setTDRStyle
from copy import deepcopy

from helpers import *
from defs import getPlot, Backgrounds, Signals, Data
import ratios

def main():
	### for data

	
	histos = ["BB","BE"]
	labels = ["dimuon_Moriond2017_BB","dimuon_Moriond2017_BE"]
	
	lambdas = [10,16,22,28,34]
	models = ["ConLL","ConLR","ConRR","DesLL","DesLR","DesRR"]
	bins = [4,7]

	massPlot = getPlot("massPlotForLimit")
	massPlotSmeared = getPlot("massPlotSmeared")
	massPlotUp = getPlot("massPlotUp")
	massPlotDown = getPlot("massPlotDown")
	massPlotWeighted = getPlot("massPlotWeighted")

	binning = [400,500,700,1100,1900,3500]


	name = "CITo2Mu_Lam10TeVConLL"
	signal = Process(getattr(Signals,name))
	
	sigHistRaw = deepcopy(signal.loadHistogramProjected(massPlot, 1))
	sigHistWeightedRaw = deepcopy(signal.loadHistogramProjected(massPlotWeighted, 1))

	sigHist = ROOT.TH1F("sigHist","sigHist",len(binning)-1,array('f',binning))
	sigHistWeighted = ROOT.TH1F("sigHistWeighted","sigHistWeighted",len(binning)-1,array('f',binning))
	sigHist.Sumw2()
	sigHistWeighted.Sumw2()	
	for index, lower in enumerate(binning):  
			if index < len(binning)-1: 
				sigHist.SetBinContent(index+1,max(0,sigHistRaw.Integral(sigHistRaw.FindBin(lower),sigHistRaw.FindBin(binning[index+1]-0.001))))
				sigHistWeighted.SetBinContent(index+1,max(0,sigHistWeightedRaw.Integral(sigHistWeightedRaw.FindBin(lower),sigHistWeightedRaw.FindBin(binning[index+1]-0.001))))
	
	
	hCanvas = TCanvas("hCanvas", "Distribution", 800,800)
	
	plotPad = ROOT.TPad("plotPad","plotPad",0,0,1,1)
	plotPad = ROOT.TPad("plotPad","plotPad",0,0.3,1,1)
	ratioPad = ROOT.TPad("ratioPad","ratioPad",0,0.,1,0.3)
	setTDRStyle()		
	plotPad.UseCurrentStyle()
	ratioPad.UseCurrentStyle()
	plotPad.Draw()	
	ratioPad.Draw()	
	plotPad.cd()
	
	plotPad.DrawFrame(400,0,3500,3500,"; m [GeV]; Events")
	
	sigHist.Draw("samehist")
	sigHistWeighted.Draw("samehist")
	sigHistWeighted.SetLineColor(kRed)
	ratioPad.cd()
	ratioGraphs =  ratios.RatioGraph(sigHist,sigHistWeighted, xMin=400, xMax=3500,title="Default/Weighted",yMin=0.8,yMax=1.2,ndivisions=10,color=ROOT.kBlack,adaptiveBinning=10000)
	ratioGraphs.draw(ROOT.gPad,True,False,True,chi2Pos=0.8)
					
	
	ROOT.gPad.RedrawAxis()
	plotPad.RedrawAxis()
	ratioPad.RedrawAxis()	
	hCanvas.Print("weighting.pdf")	
	
	name = "CITo2Mu_Lam10TeVConLL"
	signal = Process(getattr(Signals,name))
	
	sigHistRaw = deepcopy(signal.loadHistogramProjected(massPlot, 4))
	sigHistWeightedRaw = deepcopy(signal.loadHistogramProjected(massPlotWeighted, 4))

	sigHist = ROOT.TH1F("sigHist","sigHist",len(binning)-1,array('f',binning))
	sigHistWeighted = ROOT.TH1F("sigHistWeighted","sigHistWeighted",len(binning)-1,array('f',binning))
	sigHist.Sumw2()
	sigHistWeighted.Sumw2()	
	for index, lower in enumerate(binning):  
			if index < len(binning)-1: 
				sigHist.SetBinContent(index+1,max(0,sigHistRaw.Integral(sigHistRaw.FindBin(lower),sigHistRaw.FindBin(binning[index+1]-0.001))))
				sigHistWeighted.SetBinContent(index+1,max(0,sigHistWeightedRaw.Integral(sigHistWeightedRaw.FindBin(lower),sigHistWeightedRaw.FindBin(binning[index+1]-0.001))))
	
	
	hCanvas = TCanvas("hCanvas", "Distribution", 800,800)
	
	plotPad = ROOT.TPad("plotPad","plotPad",0,0,1,1)
	plotPad = ROOT.TPad("plotPad","plotPad",0,0.3,1,1)
	ratioPad = ROOT.TPad("ratioPad","ratioPad",0,0.,1,0.3)
	setTDRStyle()		
	plotPad.UseCurrentStyle()
	ratioPad.UseCurrentStyle()
	plotPad.Draw()	
	ratioPad.Draw()	
	plotPad.cd()
	
	plotPad.DrawFrame(400,0,3500,3500,"; m [GeV]; Events")
	
	sigHist.Draw("samehist")
	sigHistWeighted.Draw("samehist")
	sigHistWeighted.SetLineColor(kRed)
	ratioPad.cd()
	ratioGraphs =  ratios.RatioGraph(sigHist,sigHistWeighted, xMin=400, xMax=3500,title="Default/Weighted",yMin=0.8,yMax=1.2,ndivisions=10,color=ROOT.kBlack,adaptiveBinning=10000)
	ratioGraphs.draw(ROOT.gPad,True,False,True,chi2Pos=0.8)
					
	
	ROOT.gPad.RedrawAxis()
	plotPad.RedrawAxis()
	ratioPad.RedrawAxis()	
	hCanvas.Print("weightingBB.pdf")	
	name = "CITo2Mu_Lam10TeVConLL"
	signal = Process(getattr(Signals,name))
	
	sigHistRaw = deepcopy(signal.loadHistogramProjected(massPlot, 7))
	sigHistWeightedRaw = deepcopy(signal.loadHistogramProjected(massPlotWeighted, 7))

	sigHist = ROOT.TH1F("sigHist","sigHist",len(binning)-1,array('f',binning))
	sigHistWeighted = ROOT.TH1F("sigHistWeighted","sigHistWeighted",len(binning)-1,array('f',binning))
	sigHist.Sumw2()
	sigHistWeighted.Sumw2()
	for index, lower in enumerate(binning):  
			if index < len(binning)-1: 
				sigHist.SetBinContent(index+1,max(0,sigHistRaw.Integral(sigHistRaw.FindBin(lower),sigHistRaw.FindBin(binning[index+1]-0.001))))
				sigHistWeighted.SetBinContent(index+1,max(0,sigHistWeightedRaw.Integral(sigHistWeightedRaw.FindBin(lower),sigHistWeightedRaw.FindBin(binning[index+1]-0.001))))
	
	
	hCanvas = TCanvas("hCanvas", "Distribution", 800,800)
	
	plotPad = ROOT.TPad("plotPad","plotPad",0,0,1,1)
	plotPad = ROOT.TPad("plotPad","plotPad",0,0.3,1,1)
	ratioPad = ROOT.TPad("ratioPad","ratioPad",0,0.,1,0.3)
	setTDRStyle()		
	plotPad.UseCurrentStyle()
	ratioPad.UseCurrentStyle()
	plotPad.Draw()	
	ratioPad.Draw()	
	plotPad.cd()
	
	plotPad.DrawFrame(400,0,3500,3500,"; m [GeV]; Events")
	
	sigHist.Draw("samehist")
	sigHistWeighted.Draw("samehist")
	sigHistWeighted.SetLineColor(kRed)
	ratioPad.cd()
	ratioGraphs =  ratios.RatioGraph(sigHist,sigHistWeighted, xMin=400, xMax=3500,title="Default/Weighted",yMin=0.8,yMax=1.2,ndivisions=10,color=ROOT.kBlack,adaptiveBinning=10000)
	ratioGraphs.draw(ROOT.gPad,True,False,True,chi2Pos=0.8)
					
	
	ROOT.gPad.RedrawAxis()
	plotPad.RedrawAxis()
	ratioPad.RedrawAxis()	
	hCanvas.Print("weightingBE.pdf")	
	name = "CITo2Mu_Lam10TeVConLL"
	signal = Process(getattr(Signals,name))
	
	sigHistRaw = deepcopy(signal.loadHistogramProjected(massPlot, 10))
	sigHistWeightedRaw = deepcopy(signal.loadHistogramProjected(massPlotWeighted, 10))

	sigHist = ROOT.TH1F("sigHist","sigHist",len(binning)-1,array('f',binning))
	sigHistWeighted = ROOT.TH1F("sigHistWeighted","sigHistWeighted",len(binning)-1,array('f',binning))
	sigHist.Sumw2()
	sigHistWeighted.Sumw2()
	for index, lower in enumerate(binning):  
			if index < len(binning)-1: 
				sigHist.SetBinContent(index+1,max(0,sigHistRaw.Integral(sigHistRaw.FindBin(lower),sigHistRaw.FindBin(binning[index+1]-0.001))))
				sigHistWeighted.SetBinContent(index+1,max(0,sigHistWeightedRaw.Integral(sigHistWeightedRaw.FindBin(lower),sigHistWeightedRaw.FindBin(binning[index+1]-0.001))))
	
	
	hCanvas = TCanvas("hCanvas", "Distribution", 800,800)
	
	plotPad = ROOT.TPad("plotPad","plotPad",0,0,1,1)
	plotPad = ROOT.TPad("plotPad","plotPad",0,0.3,1,1)
	ratioPad = ROOT.TPad("ratioPad","ratioPad",0,0.,1,0.3)
	setTDRStyle()		
	plotPad.UseCurrentStyle()
	ratioPad.UseCurrentStyle()
	plotPad.Draw()	
	ratioPad.Draw()	
	plotPad.cd()
	
	plotPad.DrawFrame(400,0,3500,3500,"; m [GeV]; Events")
	
	sigHist.Draw("samehist")
	sigHistWeighted.Draw("samehist")
	sigHistWeighted.SetLineColor(kRed)
	ratioPad.cd()
	ratioGraphs =  ratios.RatioGraph(sigHist,sigHistWeighted, xMin=400, xMax=3500,title="Default/Weighted",yMin=0.8,yMax=1.2,ndivisions=10,color=ROOT.kBlack,adaptiveBinning=10000)
	ratioGraphs.draw(ROOT.gPad,True,False,True,chi2Pos=0.8)
					
	
	ROOT.gPad.RedrawAxis()
	plotPad.RedrawAxis()
	ratioPad.RedrawAxis()	
	hCanvas.Print("weightingEE.pdf")	
	
	
main()
