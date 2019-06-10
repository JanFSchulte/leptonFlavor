import argparse	
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from ROOT import TCanvas, TPad, TH1F, TH1I, THStack, TLegend, TMath, gROOT
import ratios
from setTDRStyle import setTDRStyle
gROOT.SetBatch(True)
from helpers import *
from defs import getPlot, Backgrounds, Signals, Data, Data2016, Data2018
import math
import os


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
        logbins.append(int(math.exp(math.log(m_min) + width * i)))

    return logbins

def plotDataMC(args,plot):
	

	hCanvas = TCanvas("hCanvas", "Distribution", 800,800)
	args.ratio = True
	if args.ratio:
		plotPad = ROOT.TPad("plotPad","plotPad",0,0.3,1,1)
		ratioPad = ROOT.TPad("ratioPad","ratioPad",0,0.,1,0.3)
		setTDRStyle()		
		plotPad.UseCurrentStyle()
		ratioPad.UseCurrentStyle()
		plotPad.Draw()	
		ratioPad.Draw()	
		plotPad.cd()
	else:
		plotPad = ROOT.TPad("plotPad","plotPad",0,0,1,1)
		setTDRStyle()
		plotPad.UseCurrentStyle()
		plotPad.Draw()	
		plotPad.cd()	
		
	colors = createMyColors()		

	data = Process(Data)
	data2016 = Process(Data2016)
	data2018 = Process(Data2018)
	

		
	legend = TLegend(0.55, 0.6, 0.925, 0.925)
	legend.SetFillStyle(0)
	legend.SetBorderSize(0)
	legend.SetTextFont(42)
	#~ legend.SetNColumns(2)
	legendEta = TLegend(0.15, 0.75, 0.7, 0.9)
	legendEta.SetFillStyle(0)
	legendEta.SetBorderSize(0)
	legendEta.SetTextFont(42)



	latex = ROOT.TLatex()
	latex.SetTextFont(42)
	latex.SetTextAlign(31)
	latex.SetTextSize(0.04)
	latex.SetNDC(True)
	latexCMS = ROOT.TLatex()
	latexCMS.SetTextFont(61)
	latexCMS.SetTextSize(0.06)
	latexCMS.SetNDC(True)
	latexCMSExtra = ROOT.TLatex()
	latexCMSExtra.SetTextFont(52)
	latexCMSExtra.SetTextSize(0.045)
	latexCMSExtra.SetNDC(True)	
	legendHists = []
	




	
	ROOT.gStyle.SetOptStat(0)
	
	intlumi = ROOT.TLatex()
	intlumi.SetTextAlign(12)
	intlumi.SetTextSize(0.045)
	intlumi.SetNDC(True)
	intlumi2 = ROOT.TLatex()
	intlumi2.SetTextAlign(12)
	intlumi2.SetTextSize(0.07)
	intlumi2.SetNDC(True)
	scalelabel = ROOT.TLatex()
	scalelabel.SetTextAlign(12)
	scalelabel.SetTextSize(0.03)
	scalelabel.SetNDC(True)
	metDiffLabel = ROOT.TLatex()
	metDiffLabel.SetTextAlign(12)
	metDiffLabel.SetTextSize(0.03)
	metDiffLabel.SetNDC(True)
	chi2Label = ROOT.TLatex()
	chi2Label.SetTextAlign(12)
	chi2Label.SetTextSize(0.03)
	chi2Label.SetNDC(True)
	hCanvas.SetLogy()




	
	plotPad.cd()
	plotPad.SetLogy(0)
	plotPad.SetLogy(0)
	logScale = plot.log
	
	if logScale == True:
		plotPad.SetLogy()

	datahist = data.loadHistogram(plot)	
	datahist2016 = data2016.loadHistogram(plot)	
	datahist2018 = data2018.loadHistogram(plot)	

	datahist.Scale(datahist2018.Integral(datahist2018.FindBin(60),datahist2018.FindBin(120))/datahist.Integral(datahist.FindBin(60),datahist.FindBin(120)))
	datahist2016.Scale(datahist2018.Integral(datahist2018.FindBin(60),datahist2018.FindBin(120))/datahist2016.Integral(datahist2016.FindBin(60),datahist2016.FindBin(120)))

	bng = binning()
	#~ print bng
	datahist = datahist.Rebin(len(bng) - 1, 'hist_' + uuid.uuid4().hex, array('d', bng))
	datahist2016 = datahist2016.Rebin(len(bng) - 1, 'hist_' + uuid.uuid4().hex, array('d', bng))
	datahist2018 = datahist2018.Rebin(len(bng) - 1, 'hist_' + uuid.uuid4().hex, array('d', bng))
	for i in range(0,datahist.GetNbinsX()):
		datahist.SetBinContent(i,datahist.GetBinContent(i)/datahist.GetBinWidth(i))
		datahist.SetBinError(i,datahist.GetBinError(i)/datahist.GetBinWidth(i))
		datahist2016.SetBinContent(i,datahist2016.GetBinContent(i)/datahist2016.GetBinWidth(i))
		datahist2016.SetBinError(i,datahist2016.GetBinError(i)/datahist2016.GetBinWidth(i))
		datahist2018.SetBinContent(i,datahist2018.GetBinContent(i)/datahist2018.GetBinWidth(i))
		datahist2018.SetBinError(i,datahist2018.GetBinError(i)/datahist2018.GetBinWidth(i))
	datahist2016.SetLineColor(ROOT.kBlue)
	datahist2016.SetMarkerColor(ROOT.kBlue)
	datahist.SetLineColor(ROOT.kRed)
	datahist.SetMarkerColor(ROOT.kRed)

	legend.AddEntry(datahist2016,"Data 2016 ReReco","pe")	
	legendEta.AddEntry(datahist2016,"Data 2016 ReReco","pe")	
	legend.AddEntry(datahist,"Data 2017 ReReco","pe")	
	legendEta.AddEntry(datahist,"Data 2017 ReReco","pe")	
	legend.AddEntry(datahist2018,"Data 2018 ReReco + PromptReco","pe")	
	legendEta.AddEntry(datahist2018,"Data 2018 ReReco + PromptReco","pe")	

	

	
	

	yMax = datahist.GetBinContent(datahist.GetMaximumBin())
	yMin = 0.1
	xMax = datahist.GetXaxis().GetXmax()
	xMin = datahist.GetXaxis().GetXmin()
	if plot.yMax == None:
		if logScale:
			yMax = yMax*10
		else:
			yMax = yMax*1.5
	
	else: yMax = plot.yMax
	if not plot.yMin == None:
		yMin = plot.yMin
	if not plot.xMin == None:
		xMin = plot.xMin
	if not plot.xMax == None:
		xMax = plot.xMax
	yMin = 0.00001
	plotPad.DrawFrame(xMin,yMin,xMax,yMax,"; %s ; %s" %(plot.xaxis,plot.yaxis))
	plotPad.SetLogx()
	

	datahist2016.Draw("samepe")
	datahist.Draw("samepe")
	datahist2018.Draw("samepe")




	legend.Draw()

	print datahist.Integral()

	
	latex.DrawLatex(0.95, 0.96, "36.3 - 52.7 fb^{-1} (13 TeV)")
	yLabelPos = 0.85
	cmsExtra = "Preliminary"
	#~ if not args.data:
		#~ cmsExtra = "#splitline{Private Work}{Simulation}"
		#~ yLabelPos = 0.82	
	latexCMS.DrawLatex(0.19,0.89,"CMS")
	latexCMSExtra.DrawLatex(0.19,yLabelPos,"%s"%(cmsExtra))
	


	ROOT.gPad.RedrawAxis()
	plotPad.RedrawAxis()

	if args.ratio:
		try:
			ratioPad.cd()
			ratioPad.SetLogx()
		except AttributeError:
			print "Plot fails. Look up in errs/failedPlots.txt"
			outFile =open("errs/failedPlots.txt","a")
			outFile.write('%s\n'%plot.filename%("_"+run.label+"_"+dilepton))
			outFile.close()
			plot.cuts=baseCut
			return 1
		ratioGraphs2 =  ratios.RatioGraph(datahist2016,datahist2018, xMin=xMin, xMax=xMax,title="Ratio",yMin=0.5,yMax=1.5,ndivisions=10,color=datahist2016.GetLineColor(),adaptiveBinning=100000000000)
		ratioGraphs2.draw(ROOT.gPad,True,False,True,chi2Pos=0.8)			
		ratioGraphs =  ratios.RatioGraph(datahist,datahist2018, xMin=xMin, xMax=xMax,title="Ratio",yMin=0.5,yMax=1.5,ndivisions=10,color=datahist.GetLineColor(),adaptiveBinning=100000000000)
		ratioGraphs.draw(ROOT.gPad,False,False,True,chi2Pos=0.8)



	if not os.path.exists("plots"):
		os.makedirs("plots")	
	hCanvas.Print("plots/"+plot.fileName+".pdf")

					
if __name__ == "__main__":
	
	
	parser = argparse.ArgumentParser(description='Process some integers.')
	
	parser.add_argument("-d", "--data", action="store_true", dest="data", default=False,
						  help="plot data points.")
	parser.add_argument("-m", "--mc", action="store_true", dest="mc", default=False,
						  help="plot mc backgrounds.")
	parser.add_argument("-p", "--plot", dest="plot", nargs=1, default="",
						  help="plot to plot.")
	parser.add_argument("-n", "--norm", action="store_true", dest="norm", default=False,
						  help="normalize to data.")
	parser.add_argument("-r", "--ratio", action="store_true", dest="ratio", default=False,
						  help="plot ratio plot")
	parser.add_argument("-l", "--log", action="store_true", dest="log", default=False,
						  help="plot with log scale for y axis")
	parser.add_argument("-s", "--signal", dest="signals", action="append", default=[],
						  help="signals to plot.")
	parser.add_argument("-b", "--backgrounds", dest="backgrounds", action="append", default=[],
						  help="backgrounds to plot.")


	args = parser.parse_args()
	if len(args.backgrounds) == 0:
		args.backgrounds = ["Diboson","Top","DrellYan"]

	if len(args.signals) != 0:
		args.plotSignal = True

	if args.plot == "":
		args.plot = ["massPlot"]
	
	for plot in args.plot:
		plotObject = getPlot(plot)
		plotDataMC(args,plotObject)
	
