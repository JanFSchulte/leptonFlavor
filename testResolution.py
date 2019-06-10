from ROOT import * 
from numpy import array as ar
from array import array
from setTDRStyle import setTDRStyle
from copy import deepcopy
import ratios
from helpers import *
from defs import getPlot, Backgrounds, Signals, Data

def main():

	
	resPlot1 = getPlot("resPlot1")
	resPlot2 = getPlot("resPlot2")
	resPlot3 = getPlot("resPlot3")
	resPlot4 = getPlot("resPlot4")
	resPlot5 = getPlot("resPlot5")
	resPlot6 = getPlot("resPlot6")
	resPlot7 = getPlot("resPlot7")

	drellyan = Process(getattr(Backgrounds,"DrellYan"))

	dyHist1 = deepcopy(drellyan.loadHistogram(resPlot1))
	drellyan = Process(getattr(Backgrounds,"DrellYan"))
	dyHist2 = deepcopy(drellyan.loadHistogram(resPlot2))
	drellyan = Process(getattr(Backgrounds,"DrellYan"))
	dyHist3 = deepcopy(drellyan.loadHistogram(resPlot3))
	drellyan = Process(getattr(Backgrounds,"DrellYan"))
	dyHist4 = deepcopy(drellyan.loadHistogram(resPlot4))
	drellyan = Process(getattr(Backgrounds,"DrellYan"))
	dyHist5 = deepcopy(drellyan.loadHistogram(resPlot5))
	drellyan = Process(getattr(Backgrounds,"DrellYan"))
	dyHist6 = deepcopy(drellyan.loadHistogram(resPlot6))
	drellyan = Process(getattr(Backgrounds,"DrellYan"))
	dyHist7 = deepcopy(drellyan.loadHistogram(resPlot7))
	drellyan = Process(getattr(Backgrounds,"DrellYan"))

	hCanvas = TCanvas("hCanvas", "Distribution", 800,800)

	plotPad = ROOT.TPad("plotPad","plotPad",0,0,1,1)
#	setTDRStyle()		
	plotPad.UseCurrentStyle()
	plotPad.Draw()	
	plotPad.cd()
	
	plotPad.DrawFrame(-0.5,0.000001,0.5,20000000000,"; (reco mass - gen mass)/gen mass; Events / 0.01")
	plotPad.SetLogy()
	dyHist1.Draw("samehist")
	dyHist2.Draw("samehist")
	dyHist3.Draw("samehist")
	dyHist4.Draw("samehist")
	dyHist5.Draw("samehist")
	dyHist6.Draw("samehist")
	dyHist7.Draw("samehist")
	dyHist1.SetFillColor(kWhite)
	dyHist2.SetFillColor(kWhite)
	dyHist3.SetFillColor(kWhite)
	dyHist4.SetFillColor(kWhite)
	dyHist5.SetFillColor(kWhite)
	dyHist6.SetFillColor(kWhite)
	dyHist7.SetFillColor(kWhite)
	dyHist1.SetLineColor(kRed)
	dyHist2.SetLineColor(kBlue)
	dyHist3.SetLineColor(kGreen)
	dyHist4.SetLineColor(kOrange)
	dyHist5.SetLineColor(kMagenta)
	dyHist6.SetLineColor(kCyan)
	dyHist7.SetLineColor(kBlack)
	
	legend = TLegend(0.375, 0.7, 0.925, 0.925)
	legend.SetFillStyle(0)
	legend.SetBorderSize(0)
	legend.SetTextFont(42)		
	legend.AddEntry(dyHist1,"Mass 0-250 GeV","l")	
	legend.AddEntry(dyHist2,"Mass 250-750 GeV","l")	
	legend.AddEntry(dyHist3,"Mass 750-1250 GeV","l")	
	legend.AddEntry(dyHist4,"Mass 1250-1750 GeV","l")	
	legend.AddEntry(dyHist5,"Mass 1750-2250 GeV","l")	
	legend.AddEntry(dyHist6,"Mass 2000-4000 GeV","l")	
	legend.AddEntry(dyHist7,"Mass 4000-6000 GeV","l")	
	legend.Draw()
	

	hCanvas.Print("resolution.pdf")
	
	histos = [dyHist1,dyHist2,dyHist3,dyHist4,dyHist5,dyHist6,dyHist7]
	sigmas = []
	gSystem.Load("RooCruijff_cxx.so") 	
	plotPad.SetLogy(0)

	for index, histo in enumerate(histos):
		ws = RooWorkspace("w")
		
		w = ROOT.RooWorkspace("w", ROOT.kTRUE)
		res = ROOT.RooRealVar("res","res",0.,-0.5,0.5)
		getattr(w,'import')(res)
		w.factory("weight[1.,0.,10.]")
		vars = ROOT.RooArgList(res)	
		
		dataHist1 = RooDataHist("dataHist1","dataHist1",vars, histo)
		getattr(w,'import')(dataHist1)
		w.factory("RooCruijff::cb(res, mean[0.0,-0.5,0.5], sigma[0.03,0.0,0.1], sigma[0.03,0.0,0.1], alphaL[1.43,0.,3.], alphaR[1.43,0.,3.])")
		w.pdf("cb").fitTo(w.data("dataHist1"))
		
		plotPad.DrawFrame(-0.5,0,0.5,histo.GetBinContent(histo.GetMaximumBin())*1.2,"; (reco mass - gen mass)/gen mass; Events / 0.01")
		
		frame = w.var('res').frame(RooFit.Title('mass residuals'))
		frame.GetXaxis().SetTitle('(reco mass - gen mass)/gen mass')
		frame.GetYaxis().SetTitle("Events / 0.01")
		RooAbsData.plotOn(w.data('dataHist1'), frame)
		w.pdf('cb').plotOn(frame)	
		
		frame.Draw("same")
		sigmas.append(w.var("sigma").getVal())
		hCanvas.Print("resFit%d.pdf"%index)
		
	print sigmas
	
	plotPad.DrawFrame(0,0,6000,0.1,"; mass [GeV]; mass resolution")

	graph = TGraph()
	graph.SetPoint(0,125,sigmas[0])
	graph.SetPoint(1,500,sigmas[1])
	graph.SetPoint(2,1000,sigmas[2])
	graph.SetPoint(3,1500,sigmas[3])
	graph.SetPoint(4,2000,sigmas[4])
	graph.SetPoint(5,3000,sigmas[5])
	graph.SetPoint(6,5000,sigmas[6])
	
	graph.Draw("samep")
	
	hCanvas.Print("resolutionVsMass.pdf")
		
main()
