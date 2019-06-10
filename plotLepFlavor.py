import argparse	
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from ROOT import TCanvas, TPad, TH1F, TH1I, THStack, TLegend, TMath, gROOT
import ratios
from setTDRStyle import setTDRStyle
gROOT.SetBatch(True)
from helpers import *
from defs import getPlot, Backgrounds, Backgrounds2016, Signals, Signals2016, Signals2016ADD, Data, Data2016, Data2018, path, plotList, zScale, zScale2016, zScale2018
import math
import os
from copy import copy


# Muon sys uncertainty (%) 
# as a function of mass
def getMuErr(mass, chann, norm=False):
	lumi = 0.0
	znorm = 0.0
	pileup = 0.0   # we don't use pileup 0.046
	dybkg = 0.0    # 0.07
	#pdf = 0.01*(0.433+0.003291*mass-2.159e-6*mass**2+9.044e-10*mass**3-1.807e-13*mass**4+1.51e-17*mass**5)
	pdf = 0.0
	
	# muons only next
	if chann: mutrig = 0.003
	else: mutrig = 0.007
	resolution = 0.01
	muid = 0.05

	if norm: 
		lumi = 0.0
		znorm = 0.0
		dybkg = 0
	return math.sqrt(lumi**2+znorm**2+pileup**2+dybkg**2+pdf**2+mutrig**2+resolution**2+muid**2)


# chann = True if BB
# chann = False if BE
def getElErr(mass, chann, norm=False):

	lumi = 0.0
	znorm = 0.0
	pileup = 0.046  # we don't use pileup 0.046
	dybkg = 0.0   # 0.07
	
	# poly values are in %
	#pdf = 0.01*(0.433 + 0.003291*mass - 2.159e-6*mass**2 + 9.044e-10*mass**3 - 1.807e-13*mass**4 + 1.51e-17*mass**5)
	pdf = 0.0
	
	# the following two are electrons only
	if chann: energyscale = 0.02
	else: energyscale = 0.01
	
	if chann: 
		if mass < 90: idscale = 0.01
		elif mass < 1000: idscale = 0.00002198 * mass + 0.008
		else: idscale = 0.03
	else:
		if mass < 90: idscale = 0.01
		elif mass < 300: idscale = 0.00014286 * mass - 0.00285
		else: idscale = 0.04
	
	if chann: scalefac = 0.03
	else: scalefac = 0.05
	
	if norm:
		lumi = 0.0
		znorm = 0.0
		dybkg = 0
	return math.sqrt(lumi**2+znorm**2+ pileup**2 + dybkg**2 + pdf**2 + energyscale**2 + idscale**2 + scalefac**2)


# multiply hist by 1/(Acceptance x Efficiency)
def inverseAE(hist, plotObj, year):
	# muon and electron
	# BB and BE
	if year == 2017:
		if plotObj.muon:
			if "BB" in plotObj.fileName:
				for i in range(1, hist.GetSize()-1):
					mass = hist.GetBinCenter(i)
					if mass < 600:
						ae = 2.13-0.1313*math.exp(-(mass-110.9)/20.31)-2.387*mass**(-0.03616)
					else:
						ae = 4.931-55500.0/(mass+11570.0)-0.0002108*mass
					#print mass, ae
					if mass < 120: ae = float("inf")
					hist.SetBinContent(i, hist.GetBinContent(i)*1.0/ae)
			elif "BE" in plotObj.fileName:
				for i in range(1, hist.GetSize()-1):
					mass = hist.GetBinCenter(i)
					if mass < 450:
						ae = 13.39-6.696*math.exp((mass+4855000.0)/7431000.0)-108.8*mass**(-1.138)
					else:
						ae = 0.3148+0.04447*mass**1.42*math.exp(-(mass+5108.0)/713.5)
					#print mass, ae
					if mass < 120: ae = float("inf")
					hist.SetBinContent(i, hist.GetBinContent(i)*1.0/ae)
		else: # is electron
			if "BB" in plotObj.fileName:
				for i in range(1, hist.GetSize()-1):
					mass = hist.GetBinCenter(i)
					ae = 0.5795-408.0/(mass+303.5) + 55760.0/(mass**2+98990.0)
					#print mass, ae
					hist.SetBinContent(i, hist.GetBinContent(i)*1.0/ae)
			elif "BE" in plotObj.fileName:
				for i in range(1, hist.GetSize()-1):
					mass = hist.GetBinCenter(i)
					ae = 0.01176+498.2/(mass+735.3)-100100.0/(mass**2+72990)+14190000.0/(mass**3+21600000)
					#print mass, ae
					hist.SetBinContent(i, hist.GetBinContent(i)*1.0/ae)
	elif year == 2018:
                if plotObj.muon:
                        if "BB" in plotObj.fileName:
                                for i in range(1, hist.GetSize()-1):
                                        mass = hist.GetBinCenter(i)
                                        if mass < 600:
                                                ae = 2.14-0.1286*math.exp(-(mass-110.6)/22.44)-2.366*mass**(-0.03382)
                                        else:
                                                ae = 5.18-58450.0/(mass+11570.0)-0.0002255*mass
                                        #print mass, ae
                                        if mass < 120: ae = float("inf")
                                        hist.SetBinContent(i, hist.GetBinContent(i)*1.0/ae)
                        elif "BE" in plotObj.fileName:
                                for i in range(1, hist.GetSize()-1):
                                        mass = hist.GetBinCenter(i)
                                        if mass < 450:
                                                ae = 13.4-6.693*math.exp((mass+4852000.0)/7437000.0)-81.43*mass**(-1.068)
                                        else:
                                                ae = 0.3154+0.04561*mass**1.362*math.exp(-(mass+4927.0)/727.5)
                                        #print mass, ae
                                        if mass < 120: ae = float("inf")
                                        hist.SetBinContent(i, hist.GetBinContent(i)*1.0/ae)
                else: # is electron
                        if "BB" in plotObj.fileName:
                                for i in range(1, hist.GetSize()-1):
                                        mass = hist.GetBinCenter(i)
                                        ae = 0.5947-440.1/(mass+393) + 47630.0/(mass**2+108000)
                                        #print mass, ae
                                        hist.SetBinContent(i, hist.GetBinContent(i)*1.0/ae)
                        elif "BE" in plotObj.fileName:
                                for i in range(1, hist.GetSize()-1):
                                        mass = hist.GetBinCenter(i)
                                        ae = 0.01718+468.9/(mass+575.6)-113300.0/(mass**2+82800)+13740000.0/(mass**3+23380000)
                                        #print mass, ae
                                        hist.SetBinContent(i, hist.GetBinContent(i)*1.0/ae)

		

def plotDataMC(args,plot_mu,plot_el):
	

	hCanvas = TCanvas("hCanvas", "Distribution", 800,800)
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
		
	# Data load processes
	colors = createMyColors()		
	if args.use2016:
		data = Process(Data2016, normalized=True)
	elif args.use2018:
		data = Process(Data2018, normalized=True)
	else:	
		data = Process(Data, normalized=True)
	
	eventCounts_mu = totalNumberOfGeneratedEvents(path,plot_mu.muon)	
	eventCounts_el = totalNumberOfGeneratedEvents(path,plot_el.muon)
	negWeights_mu = negWeightFractions(path,plot_mu.muon)
	negWeights_el = negWeightFractions(path,plot_el.muon)
	#print negWeights

	# Background load processes	
	backgrounds = copy(args.backgrounds)
	'''if plot_mu.useJets:
		if "Wjets" in backgrounds:
			backgrounds.remove("Wjets")
		backgrounds.insert(0,"Jets")'''
	processes_mu = []
	processes_el = []
	for background in backgrounds:
		if args.use2016:
			if background == "Jets":
				processes_mu.append(Process(getattr(Backgrounds2016,background),eventCounts_mu,negWeights_mu,normalized=True))
				processes_el.append(Process(getattr(Backgrounds2016,background),eventCounts_el,negWeights_el,normalized=True))
			else:	
				processes_mu.append(Process(getattr(Backgrounds2016,background),eventCounts_mu,negWeights_mu))
				processes_el.append(Process(getattr(Backgrounds2016,background),eventCounts_el,negWeights_el))
		else:
			if background == "Jets":
				processes_mu.append(Process(getattr(Backgrounds,background),eventCounts_mu,negWeights_mu,normalized=True))
				processes_el.append(Process(getattr(Backgrounds,background),eventCounts_el,negWeights_el,normalized=True))
			else:	
				processes_mu.append(Process(getattr(Backgrounds,background),eventCounts_mu,negWeights_mu))
				processes_el.append(Process(getattr(Backgrounds,background),eventCounts_el,negWeights_el))

	
	'''# Signal load processes
	signals = []
	for signal in args.signals:
		if args.use2016:
			if args.ADD: signals.append(Process(getattr(Signals2016ADD, signal), eventCounts, negWeights))
			else: signals.append(Process(getattr(Signals2016,signal),eventCounts,negWeights))
		else:	
			if args.ADD: signals.append(Process(getattr(SignalsADD, signal), eventCounts, negWeights))
			else: signals.append(Process(getattr(Signals,signal),eventCounts,negWeights))
	'''	
	legend = TLegend(0.55, 0.75, 0.925, 0.925)
	legend.SetFillStyle(0)
	legend.SetBorderSize(0)
	legend.SetTextFont(42)
	
	'''legendEta = TLegend(0.35, 0.55, 0.9, 0.9)
	legendEta.SetFillStyle(0)
	legendEta.SetBorderSize(0)
	legendEta.SetTextFont(42)
	legendEta.SetNColumns(2)
	'''

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
	
	# Modify legend information
	legendHistData = ROOT.TH1F()
	if args.data:	
		legend.AddEntry(legendHistData,"Data","pe")	
		legendEta.AddEntry(legendHistData,"Data","pe")	
	
	for process in reversed(processes_mu):
		if not plot_mu.muon and "#mu^{+}#mu^{-}" in process.label:
			process.label = process.label.replace("#mu^{+}#mu^{-}","e^{+}e^{-}")
		process.theColor = ROOT.kBlue
		process.theLineColor = ROOT.kBlue
		temphist = ROOT.TH1F()
		temphist.SetFillColor(process.theColor)
		legendHists.append(temphist.Clone)
		legend.AddEntry(temphist,process.label,"f")
		#legendEta.AddEntry(temphist,process.label,"f")
	
        for process in reversed(processes_el):
                if not plot_el.muon and "#mu^{+}#mu^{-}" in process.label:
                        process.label = process.label.replace("#mu^{+}#mu^{-}","e^{+}e^{-}")
		process.theColor = ROOT.kRed
		process.theLineColor = ROOT.kRed
                temphist = ROOT.TH1F()
                temphist.SetFillColor(process.theColor)
                legendHists.append(temphist.Clone)
                legend.AddEntry(temphist,process.label,"f")
                #legendEta.AddEntry(temphist,process.label,"f")

	'''if args.signals !=0:
		processesWithSignal = []
		for process in processes:
			processesWithSignal.append(process)
		for Signal in signals:
			processesWithSignal.append(Signal)
			temphist = ROOT.TH1F()
			temphist.SetFillColor(Signal.theColor)
			temphist.SetLineColor(Signal.theLineColor)
			legendHists.append(temphist.Clone)		
			legend.AddEntry(temphist,Signal.label,"l")
			legendEta.AddEntry(temphist,Signal.label,"l")
	'''

	# Modify plot pad information	
	nEvents=-1

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


	# Luminosity information	
	plotPad.cd()
	plotPad.SetLogy(0)
	logScale = plot_mu.log
	
	if logScale == True:
		plotPad.SetLogy()

	if args.use2016:	
		lumi_el = 35.9*1000
		lumi_mu = 36.3*1000
	elif args.use2018:	
		lumi_el = 59.97*1000
		lumi_mu = 61.608*1000
	else:
		lumi_el = 41.529*1000
		lumi_mu = 42.135*1000
	if args.use2016:		
		zScaleFac_mu = zScale2016["muons"]
		zScaleFac_el = zScale2016["electrons"]
	elif args.use2018:		
		zScaleFac_mu = zScale2018["muons"]
		zScaleFac_el = zScale2018["electrons"]
	else:
		zScaleFac_mu = zScale["muons"]
		zScaleFac_el = zScale["electrons"]
	
			
	# Data and background loading	
	'''if plot.plot2D:	
		datahist = data.loadHistogramProjected(plot,lumi,zScaleFac)	
		
		stack = TheStack2D(processes,lumi,plot,zScaleFac)
	else:
		datahist = data.loadHistogram(plot,lumi,zScaleFac)	
		
		stack = TheStack(processes,lumi,plot,zScaleFac)
	'''
	lumi_mu = 1.0 * 1000
	lumi_el = 1.0 * 1000
	datahist = data.loadHistogram(plot_mu,lumi_mu,zScaleFac_mu)
	stackmu = TheStack(processes_mu,lumi_mu,plot_mu,zScaleFac_mu)
	stackel = TheStack(processes_el,lumi_el,plot_el,zScaleFac_el)
	# call hist in stack by: for h in stackmu.theStack.GetHists()
	
	#muheight = stackmu.theHistogram.FindBin(90)
	#print "Z height in mu: %.3f, %.3f"%(stackmu.theHistogram.GetBinContent(muheight), stackel.theHistogram.GetBinContent(muheight))
	if args.znorm:
		muheight = stackmu.theHistogram.FindBin(90)
		print "Z height of mu: %d +- %d"%(stackmu.theHistogram.GetBinCenter(muheight), stackmu.theHistogram.GetBinWidth(muheight))
		print "Z height of mu: %d"%(stackmu.theHistogram.GetBinContent(muheight))
		elheight = stackel.theHistogram.FindBin(90)
		print "Z height of el: %d +- %d"%(stackel.theHistogram.GetBinCenter(elheight), stackel.theHistogram.GetBinWidth(elheight))
		print "Z height of el: %d"%(stackel.theHistogram.GetBinContent(elheight))
		znum = stackmu.theHistogram.GetBinContent(muheight)
		print znum
		for h in stackmu.theStack.GetHists(): h.Scale(1./znum)
		for h in stackel.theStack.GetHists(): h.Scale(1./znum)
		stackmu.theHistogram.Scale(1./znum)
		stackel.theHistogram.Scale(1./znum)
	
	if args.ae:
		year = 2017
		if args.use2018: year = 2018
		for h in stackmu.theStack.GetHists(): inverseAE(h, plot_mu, year)
		for h in stackel.theStack.GetHists(): inverseAE(h, plot_el, year)
		inverseAE(stackmu.theHistogram, plot_mu, year)
		inverseAE(stackel.theHistogram, plot_el, year)
	
	if args.data:
		yMax = datahist.GetBinContent(datahist.GetMaximumBin())
		if "Mass" in plot.fileName:
			yMin = 0.00001
		else:
			yMin = 0.01
		xMax = datahist.GetXaxis().GetXmax()
		xMin = datahist.GetXaxis().GetXmin()
	else:	
		yMax = stackmu.theHistogram.GetBinContent(datahist.GetMaximumBin())
		yMin = 0.01
		xMax = stackmu.theHistogram.GetXaxis().GetXmax()
		xMin = stackmu.theHistogram.GetXaxis().GetXmin()	
	if plot_mu.yMax == None:
		if logScale:
			yMax = yMax*10000
		else:
			yMax = yMax*1.5
	else: yMax = plot_mu.yMax
	
	if "Mass" in plot_mu.fileName:
		yMax = 20000000	
	
	if not plot_mu.yMin == None:
		yMin = plot_mu.yMin
	if not plot_mu.xMin == None:
		xMin = plot_mu.xMin
	if not plot_mu.xMax == None:
		xMax = plot_mu.xMax
	#if args.ADD and args.use2016: 
	#	xMin = 1700
	#	xMax = 4000
	#	yMax = 1.0
	#if "CosThetaStarBBM1800" in plot.fileName:
	#	yMax = 3
	#print xMin, xMax, yMin, yMax
	'''xMin = 70
	xMax = 4000
	yMin = 0.00001
	yMax = 20000000'''
	yMin = 0.00001 / 40
	yMax = 200000000.0 / 40
	if args.ae: xMin = 200
	if args.znorm: 
		yMin /= 300
		yMax /= 300
	plotPad.DrawFrame(xMin,yMin,xMax,yMax,"; %s ; %s" %("m(l^{+}l^{-}) [GeV]","d#sigma(pp#rightarrow ll)"))
	
	
	drawStack_mu = stackmu
	drawStack_el = stackel
 	#~ print datahist.Integral(datahist.FindBin(60),datahist.FindBin(120))/drawStack.theHistogram.Integral(drawStack.theHistogram.FindBin(60),drawStack.theHistogram.FindBin(120))
 	#~ low = 900
 	#~ high = 1300
 	#~ print datahist.Integral(datahist.FindBin(low),datahist.FindBin(high))
 	#~ print drawStack.theHistogram.Integral(datahist.FindBin(low),datahist.FindBin(high))

	
	# Draw background from stack
	drawStack_mu.theStack.Draw("same hist")
	drawStack_el.theStack.Draw("same hist")


	'''# Draw signal information
	if len(args.signals) != 0:
		signalhists = []
		for Signal in signals:
			if plot.plot2D: # plot collins-soper angle
				signalhist = Signal.loadHistogramProjected(plot,lumi, zScaleFac)
				signalhist.SetLineWidth(2)
				signalBackgrounds = deepcopy(backgrounds)
				signalBackgrounds.remove("DrellYan")
				signalProcesses = []
				for background in signalBackgrounds:
					if background == "Jets":
						signalProcesses.append(Process(getattr(Backgrounds,background),eventCounts,negWeights,normalized=True))
					else:	
						signalProcesses.append(Process(getattr(Backgrounds,background),eventCounts,negWeights))
				signalStack = TheStack2D(signalProcesses,lumi,plot, zScaleFac)
				signalhist.Add(signalStack.theHistogram)
				signalhist.SetMinimum(0.1)
				signalhist.Draw("samehist")
				signalhists.append(signalhist)	
			else:
				signalhist = Signal.loadHistogram(plot,lumi,zScaleFac)
				signalhist.SetLineWidth(2)
				signalBackgrounds = deepcopy(backgrounds)
				signalBackgrounds.remove("DrellYan") # signalBackgrounds = ["Jets", "Other"]
				signalProcesses = []
				for background in signalBackgrounds:
					if background == "Jets":
						signalProcesses.append(Process(getattr(Backgrounds,background),eventCounts,negWeights,normalized=True))
					else:	
						signalProcesses.append(Process(getattr(Backgrounds,background),eventCounts,negWeights))
				signalStack = TheStack(signalProcesses,lumi,plot,zScaleFac)
				signalhist.Add(signalStack.theHistogram)
				signalhist.SetMinimum(0.0001)
				signalhist.Draw("samehist")
				signalhists.append(signalhist)	
	'''
	# Draw data
	datahist.SetMinimum(0.0001)
	#if args.data:
	#	datahist.Draw("samep")	

	# Draw legend
	if "Eta" in plot_mu.fileName or "CosTheta" in plot_mu.fileName:
		legendEta.Draw()
	else:
		legend.Draw()

	plotPad.SetLogx(plot_mu.logX)
	
	latex.DrawLatex(0.95, 0.96, "%.3f fb^{-1} (13 TeV, #mu#mu), %.3f fb^{-1} (13 TeV, ee)"%(lumi_mu*0.001, lumi_el*0.001))
	yLabelPos = 0.85
	cmsExtra = "Private Work"
	if not args.data:
		cmsExtra = "#splitline{Private Work}{Simulation}"
		yLabelPos = 0.82	
	latexCMS.DrawLatex(0.19,0.89,"CMS")
	latexCMSExtra.DrawLatex(0.19,yLabelPos,"%s"%(cmsExtra))
	#~ print datahist.Integral()
	if args.ratio:
		try:
			ratioPad.cd()
			ratioPad.SetLogx(plot_mu.logX)
		except AttributeError:
			print ("Plot fails. Look up in errs/failedPlots.txt")
			outFile =open("errs/failedPlots.txt","a")
			outFile.write('%s\n'%plot_mu.filename%("_"+run.label+"_"+dilepton))
			outFile.close()
			plot_mu.cuts=baseCut
			return 1
		#ratioGraphs = ratios.RatioGraph(drawStack_mu.theStack.GetHists()[0], drawStack_el.theStack.GetHists()[0], xMin=xMin, xMax=xMax,title="R_{#mu#mu/ee}",yMin=0.0,yMax=2,ndivisions=10,color=ROOT.kBlack,adaptiveBinning=10000)
		#setRatioError(ratioGraphs, drawStack_mu.theHistogram, drawStack_el.theHistogram)
		#ratioGraphs.draw(ROOT.gPad,True,False,True,chi2Pos=0.8)
		hhmu = drawStack_mu.theHistogram
		hhel = drawStack_el.theHistogram
		ratioGraphs = ROOT.TGraphAsymmErrors(hhmu.GetSize()-2)
		chann = True if "BB" in plot_mu.fileName else False
		for i in range(1, hhmu.GetSize()-1):
			xval = hhmu.GetBinCenter(i)
			xerr = hhmu.GetBinWidth(i)/2
			if hhel.GetBinContent(i) == 0: continue
			if hhel.GetBinContent(i) < 0 or hhmu.GetBinContent(i) < 0: continue
			yval = hhmu.GetBinContent(i)*1.0/hhel.GetBinContent(i)
			#if yval > 10: continue
			yerr = yval * math.sqrt(getMuErr(xval, chann, args.znorm)**2 + getElErr(xval, chann, args.znorm)**2)
			ratioGraphs.SetPoint(i, xval, yval)
			ratioGraphs.SetPointError(i, xerr, xerr, yerr, yerr)
		nBinsX = 20
                nBinsY = 10
                hAxis = ROOT.TH2F("hAxis", "", nBinsX, xMin, xMax, nBinsY, 0.0, 4.0)
                hAxis.Draw("AXIS")

                hAxis.GetYaxis().SetNdivisions(408)
                hAxis.SetTitleOffset(0.4, "Y")
                hAxis.SetTitleSize(0.15, "Y")
                hAxis.SetYTitle("R_{#mu#mu/ee}")
                hAxis.GetXaxis().SetLabelSize(0.0)
                hAxis.GetYaxis().SetLabelSize(0.15)
		hAxis.SetTitleSize(0.15, "Y")
                #binMerging = [-1]
		
                oneLine = ROOT.TLine(xMin, 1.0, xMax, 1.0)
                oneLine.SetLineStyle(2)
                oneLine.Draw()
		
		ratioGraphs.SetFillColor(6)
		ratioGraphs.SetFillStyle(3002)	
		#ratioGraphs.Draw("SAMEpZ4")
		ratioGraphs.Draw("same p")
		ratioPad.Update()
					

	ROOT.gPad.RedrawAxis()
	plotPad.RedrawAxis()
	if args.ratio:
		ratioPad.RedrawAxis()
	if not os.path.exists("lepFlavor"):
		os.makedirs("lepFlavor")	

	if args.use2016: year = "2016"
	elif args.use2018: year = "2018"
	else: year = "2017"

	if args.ae: year += "_inverseAE"
	if args.znorm: year += "_znorm"
	
	hCanvas.Print("lepFlavor/%s_%s_other.pdf"%(plot_mu.fileName, year))


					
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
	parser.add_argument("-2016", "--2016", action="store_true", dest="use2016", default=False,
						  help="use 2016 data and MC.")
	parser.add_argument("-2018", "--2018", action="store_true", dest="use2018", default=False,
						  help="use 2018 data with 2017 MC.")
	parser.add_argument("-r", "--ratio", action="store_true", dest="ratio", default=False,
						  help="plot ratio plot")
	parser.add_argument("-l", "--log", action="store_true", dest="log", default=False,
						  help="plot with log scale for y axis")
	parser.add_argument("-s", "--signal", dest="signals", action="append", default=[],
						  help="signals to plot.")
	parser.add_argument("-b", "--backgrounds", dest="backgrounds", action="append", default=[],
						  help="backgrounds to plot.")
	parser.add_argument("-a", "--ADD", action="store_true", dest="ADD", default=False, help="plot add signals")
	parser.add_argument("--ae", action="store_true", dest="ae", default=False,help="times inverse Acceptance x Efficiency")
	parser.add_argument("--znorm", action="store_true", dest="znorm", default=False, help="normalize to z peak")

	args = parser.parse_args()
	if len(args.backgrounds) == 0:
		#args.backgrounds = ["Wjets","Other","DrellYan"]
		#~ args.backgrounds = ["Diboson","DrellYan"]
		args.backgrounds = ["Other"]

	if len(args.signals) != 0:
		args.plotSignal = True

	'''if args.plot == "":
		args.plot = plotList
	'''
	muplots = ["massPlotBB", "massPlotBE"]
	elplots = ["massPlotEleBB", "massPlotEleBE"]
	signals = args.signals
	for i in range(len(muplots)):
		args.signals = signals
		plot_mu = getPlot(muplots[i])
		plot_el = getPlot(elplots[i])
		#plot_mu.logX = False
		#plot_el.logX = False
		'''if len(args.signals) > 0:
			#~ if ("To2E" in args.signals[0] and plotObject.muon) or ("To2Mu" in args.signals[0] and not plotObject.muon):
			args.signals = []
			if plotObject.muon:
				for signal in signals:
					if args.ADD: args.signals.append("ADDGravTo2Mu_"+signal)
					else: args.signals.append("CITo2Mu_"+signal)
			else:
				for signal in signals:
					if args.ADD: args.signals.append("ADDGravTo2E_"+signal)
					else: args.signals.append("CITo2E_"+signal)
		'''
		#~ print args.plotSignal	
		plotDataMC(args,plot_mu,plot_el)
	
