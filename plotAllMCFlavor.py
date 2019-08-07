import argparse	
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from ROOT import TCanvas, TPad, TH1F, TH1I, THStack, TLegend, TMath, gROOT, TGaxis
import ratios
from setTDRStyle import setTDRStyle
gROOT.SetBatch(True)
from helpers import *
from defs import getPlot, Backgrounds, Backgrounds2016, Backgrounds2018, Signals, Signals2016, Signals2016ADD, Data, Data2016, Data2018, path, plotList, zScale, zScale2016, zScale2018
import math
import os
from copy import copy
import numpy as np
import root_numpy

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

def getErrors(default, others):
	dfarr=root_numpy.hist2array(default)
	errs=np.zeros(len(dfarr))
	for other in others:
		if type(other)==list:
			err1=root_numpy.hist2array(other[0])-dfarr
			err1=err1**2
			err2=root_numpy.hist2array(other[1])-dfarr
			err2=err2**2
			err=np.maximum(err1,err2)
			errs+=err
		else:           
			err=root_numpy.hist2array(other)-dfarr
			err=err**2
			errs+=err
	return errs                           

# multiply hist by 1/(Acceptance x Efficiency)
def inverseAE(hist, plotObj, year):
		# muon and electron
		# BB and BE
	if year == 2016:
		if plotObj.muon:
			if "BB" in plotObj.fileName:
				for i in range(1, hist.GetSize()-1):
					mass = hist.GetBinCenter(i)
					if mass < 600:
						ae = 2.129-0.1268*math.exp(-(mass-119.2)/22.35)-2.386*mass**(-0.03619)
					else:
						ae = 2.891-2.291e+04/(mass+8294.)-0.0001247*mass
					#print mass, ae
					if mass < 120: ae = float("inf")
					hist.SetBinContent(i, hist.GetBinContent(i)*1.0/ae)
			elif "BE" in plotObj.fileName:
				for i in range(1, hist.GetSize()-1):
					mass = hist.GetBinCenter(i)
					if mass < 450:
									ae = 13.56-6.672*math.exp((mass+4.879e+06)/7.233e+06)-826*mass**(-1.567)
					else:
									ae =  0.2529+0.06511*mass**0.8755*math.exp(-(mass+4601.)/1147)
					#print mass, ae
					if mass < 120: ae = float("inf")
					hist.SetBinContent(i, hist.GetBinContent(i)*1.0/ae)
		else: # is electron
			if "BB" in plotObj.fileName:
				for i in range(1, hist.GetSize()-1):
					mass = hist.GetBinCenter(i)
					ae = 0.6386-497.7/(mass+348.7) + 69570.0/(mass**2+115400.0)
					#print mass, ae
					hist.SetBinContent(i, hist.GetBinContent(i)*1.0/ae)
			elif "BE" in plotObj.fileName:
				for i in range(1, hist.GetSize()-1):
					mass = hist.GetBinCenter(i)
					ae = -0.03377+735.1/(mass+1318)-88890.0/(mass**2+75720)+14240000.0/(mass**3+23420000)
					#print mass, ae
					hist.SetBinContent(i, hist.GetBinContent(i)*1.0/ae)

	elif year == 2017:
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
								ae = 0.585-404.6/(mass+279.5) + 56180.0/(mass**2+91430.0)
								#print mass, ae
								hist.SetBinContent(i, hist.GetBinContent(i)*1.0/ae)
				elif "BE" in plotObj.fileName:
						for i in range(1, hist.GetSize()-1):
								mass = hist.GetBinCenter(i)
								ae = 0.02066+429.7/(mass+729)-90960.0/(mass**2+71900)+13780000.0/(mass**3+22050000)
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
									ae = 0.576-417.7/(mass+381.8) + 46070.0/(mass**2+107200)
									#print mass, ae
									hist.SetBinContent(i, hist.GetBinContent(i)*1.0/ae)
					elif "BE" in plotObj.fileName:
							for i in range(1, hist.GetSize()-1):
									mass = hist.GetBinCenter(i)
									ae = 0.01443+475.7/(mass+639.1)-105600.0/(mass**2+82810)+12890000.0/(mass**3+23170000)
									#print mass, ae
									hist.SetBinContent(i, hist.GetBinContent(i)*1.0/ae)
																				
def Stacks(processes,lumi,plot,zScale):
	stacks=[]
	for i in range(3):
		stacks.append(TheStack(processes[i],lumi[i],plot,zScale[i]))
	return stacks
def Addhist(histlist):
	tempHist=histlist[0]
	for i in range(1,3):
		tempHist.Add(histlist[i])
	return tempHist         
def Addstack(Stacklist):
	tempStack=Stacklist[0]
	for i in range(1,3):
		tempStack.Add(Stacklist[i])
	return tempStack                                                                                                                                                                          
def plotDataMC(args,plot_mu,plot_el):
	

	hCanvas = TCanvas("hCanvas", "Distribution", 800,800)
	if args.ratio:
		plotPad = ROOT.TPad("plotPad","plotPad",0,0.5,1,1)
		ratioPad = ROOT.TPad("ratioPad","ratioPad",0,0.,1,0.5)
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
		data_mu = Process(Data2016, normalized=True)
		data_el = Process(Data2016, normalized=True)
	elif args.use2018:
		data_mu = Process(Data2018, normalized=True)
		data_el = Process(Data2018, normalized=True)
	elif args.useall:
		data_all=[Process(Data2016, normalized=True),Process(Data, normalized=True),Process(Data2018, normalized=True)]              
	else:	
		data_mu = Process(Data, normalized=True)
		data_el = Process(Data, normalized=True)
	
	eventCounts_mu = totalNumberOfGeneratedEvents(path,plot_mu["default"].muon)	
	eventCounts_el = totalNumberOfGeneratedEvents(path,plot_el["default"].muon)
	negWeights_mu = negWeightFractions(path,plot_mu["default"].muon)
	negWeights_el = negWeightFractions(path,plot_el["default"].muon)

	# Background load processes	
	backgrounds = copy(args.backgrounds)
	if plot_mu["default"].useJets:
		if "Wjets" in backgrounds:
			backgrounds.remove("Wjets")
		backgrounds.insert(0,"Jets")
	processes_mu = []
	processes_el = []
	processes_mu2016=[]
	processes_mu2017=[]
	processes_mu2018=[]
	processes_el2016=[]
	processes_el2017=[]
	processes_el2018=[]
	for background in backgrounds:
		if args.use2016:
			if background == "Jets":
				processes_mu.append(Process(getattr(Backgrounds2016,background),eventCounts_mu,negWeights_mu,normalized=True))
				processes_el.append(Process(getattr(Backgrounds2016,background),eventCounts_el,negWeights_el,normalized=True))
			else:	
				processes_mu.append(Process(getattr(Backgrounds2016,background),eventCounts_mu,negWeights_mu))
				processes_el.append(Process(getattr(Backgrounds2016,background),eventCounts_el,negWeights_el))
		elif args.use2018:
				if background == "Jets":
						processes_mu.append(Process(getattr(Backgrounds2018,background),eventCounts_mu,negWeights_mu,normalized=True))
						processes_el.append(Process(getattr(Backgrounds2018,background),eventCounts_el,negWeights_el,normalized=True))
				else:
						processes_mu.append(Process(getattr(Backgrounds2018,background),eventCounts_mu,negWeights_mu))
						processes_el.append(Process(getattr(Backgrounds2018,background),eventCounts_el,negWeights_el))
		elif args.useall:
				if background == "Jets":
						processes_mu2016.append(Process(getattr(Backgrounds2016,background),eventCounts_mu,negWeights_mu,normalized=True))
						processes_el2016.append(Process(getattr(Backgrounds2016,background),eventCounts_el,negWeights_el,normalized=True))
						processes_mu2017.append(Process(getattr(Backgrounds,background),eventCounts_mu,negWeights_mu,normalized=True))
						processes_el2017.append(Process(getattr(Backgrounds,background),eventCounts_el,negWeights_el,normalized=True))
						processes_mu2018.append(Process(getattr(Backgrounds2018,background),eventCounts_mu,negWeights_mu,normalized=True))
						processes_el2018.append(Process(getattr(Backgrounds2018,background),eventCounts_el,negWeights_el,normalized=True))
						processes_mu=[processes_mu2016,processes_mu2017,processes_mu2018]
						processes_el=[processes_mu2016,processes_mu2017,processes_mu2018]
				else:
						processes_mu2016.append(Process(getattr(Backgrounds2016,background),eventCounts_mu,negWeights_mu))
						processes_el2016.append(Process(getattr(Backgrounds2016,background),eventCounts_el,negWeights_el))
						processes_mu2017.append(Process(getattr(Backgrounds,background),eventCounts_mu,negWeights_mu))
						processes_el2017.append(Process(getattr(Backgrounds,background),eventCounts_el,negWeights_el))
						processes_mu2018.append(Process(getattr(Backgrounds2018,background),eventCounts_mu,negWeights_mu))
						processes_el2018.append(Process(getattr(Backgrounds2018,background),eventCounts_el,negWeights_el))
						processes_mu=[processes_mu2016,processes_mu2017,processes_mu2018]
						processes_el=[processes_mu2016,processes_mu2017,processes_mu2018]
		else:             
				if background == "Jets":
						processes_mu.append(Process(getattr(Backgrounds,background),eventCounts_mu,negWeights_mu,normalized=True))
						processes_el.append(Process(getattr(Backgrounds,background),eventCounts_el,negWeights_el,normalized=True))
				else:
						processes_mu.append(Process(getattr(Backgrounds,background),eventCounts_mu,negWeights_mu))
						processes_el.append(Process(getattr(Backgrounds,background),eventCounts_el,negWeights_el))

	
	legend = TLegend(0.55, 0.75, 0.925, 0.925)
	legend.SetFillStyle(0)
	legend.SetBorderSize(0)
	legend.SetTextFont(42)
	

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
	legendHistData_mu = ROOT.TH1F()
	legendHistData_el = ROOT.TH1F()
	dy_mu = ROOT.TH1F()
	dy_el = ROOT.TH1F()
	if args.data:	
		legendHistData_mu.SetMarkerColor(ROOT.kViolet)
		legendHistData_el.SetMarkerColor(ROOT.kOrange)
		dy_mu.SetLineColor(ROOT.kBlue-3)
		dy_el.SetLineColor(ROOT.kRed-3)
		legend.AddEntry(legendHistData_mu,"Data #rightarrow #mu^{+}#mu^{-}","pe")
		legend.AddEntry(legendHistData_el,"Data #rightarrow e^{+}e^{-}", "pe")
		legend.AddEntry(dy_mu, "MC Inclusive #rightarrow #mu^{+}#mu^{-}", "l")
		legend.AddEntry(dy_el, "MC Inclusive #rightarrow e^{+}e^{-}", "l")
		 	
	if args.useall:
		for i in range(3):
			for process in reversed(processes_mu[i]):
				if not plot_mu["default"].muon and "#mu^{+}#mu^{-}" in process.label:
					process.label = process.label.replace("#mu^{+}#mu^{-}","e^{+}e^{-}")
				process.theColor = ROOT.kBlue
				process.theLineColor = ROOT.kBlue
				temphist = ROOT.TH1F()
				temphist.SetFillColor(process.theColor)

			for process in reversed(processes_el[i]):
				if not plot_el["default"].muon and "#mu^{+}#mu^{-}" in process.label:
					process.label = process.label.replace("#mu^{+}#mu^{-}","e^{+}e^{-}")
				process.theColor = ROOT.kRed
				process.theLineColor = ROOT.kRed
				temphist = ROOT.TH1F()
				temphist.SetFillColor(process.theColor)
	else:
		for process in reversed(processes_mu):
			if not plot_mu["default"].muon and "#mu^{+}#mu^{-}" in process.label:
				process.label = process.label.replace("#mu^{+}#mu^{-}","e^{+}e^{-}")
			process.theColor = ROOT.kBlue
			process.theLineColor = ROOT.kBlue
			temphist = ROOT.TH1F()
			temphist.SetFillColor(process.theColor)
	
		for process in reversed(processes_el):
			if not plot_el["default"].muon and "#mu^{+}#mu^{-}" in process.label:
				process.label = process.label.replace("#mu^{+}#mu^{-}","e^{+}e^{-}")
			process.theColor = ROOT.kRed
			process.theLineColor = ROOT.kRed
			temphist = ROOT.TH1F()
			temphist.SetFillColor(process.theColor)

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
	logScale = plot_mu["default"].log
	
	if logScale == True:
		plotPad.SetLogy()

	if args.use2016:	
		lumi_el = 35.9*1000
		lumi_mu = 36.3*1000
	elif args.use2018:	
		lumi_el = 59.97*1000
		lumi_mu = 61.608*1000
	elif args.useall:
		lumi_el = [35.9*1000,41.529*1000,59.97*1000]
		lumi_mu = [36.3*1000,42.135*1000,61.608*1000]
	else:
		lumi_el = 41.529*1000
		lumi_mu = 42.135*1000
	if args.use2016:		
		zScaleFac_mu = zScale2016["muons"]
		zScaleFac_el = zScale2016["electrons"]
	elif args.use2018:		
		zScaleFac_mu = zScale2018["muons"]
		zScaleFac_el = zScale2018["electrons"]
	elif args.useall:
		zScaleFac_mu = [zScale2016["muons"],zScale["muons"],zScale2018["muons"]]
		zScaleFac_el = [zScale2016["electrons"],zScale["electrons"],zScale2018["electrons"]]
	else:	
		zScaleFac_mu = zScale["muons"]
		zScaleFac_el = zScale["electrons"]
	
			
	# Data and background loading
	if args.useall:
		datamu=[]
		datael=[]
		for i in range(3):
			datamu.append(data_all[i].loadHistogram(plot_mu["default"],lumi_mu[i],zScaleFac_mu[i]))
			datael.append(data_all[i].loadHistogram(plot_el["default"],lumi_el[i],zScaleFac_el[i]))
		stackmu = Stacks(processes_mu,lumi_mu,plot_mu["default"],zScaleFac_mu)
		mu_scaleup=Stacks(processes_mu,lumi_mu,plot_mu["scale_up"],zScaleFac_mu)
		mu_scaledown=Stacks(processes_mu,lumi_mu,plot_mu["scale_down"],zScaleFac_mu)
		mu_ID=Stacks(processes_mu,lumi_mu,plot_mu["ID"],zScaleFac_mu)
		mu_reso=Stacks(processes_mu,lumi_mu,plot_mu["reso"],zScaleFac_mu)
		stackel = Stacks(processes_el,lumi_el,plot_el["default"],zScaleFac_el)
		el_scaleup=Stacks(processes_el,lumi_el,plot_el["scale_up"],zScaleFac_el)
		el_scaledown=Stacks(processes_el,lumi_el,plot_el["scale_down"],zScaleFac_el)
		el_PUup=Stacks(processes_el,lumi_el,plot_el["PU_up"],zScaleFac_el)
		el_PUdown=Stacks(processes_el,lumi_el,plot_el["PU_down"],zScaleFac_el)
	else:	
		datamu = data_mu.loadHistogram(plot_mu["default"],lumi_mu,zScaleFac_mu)
		datael = data_el.loadHistogram(plot_el["default"],lumi_el,zScaleFac_el)
		stackmu = TheStack(processes_mu,lumi_mu,plot_mu["default"],zScaleFac_mu)
		mu_scaleup=TheStack(processes_mu,lumi_mu,plot_mu["scale_up"],zScaleFac_mu)
		mu_scaledown=TheStack(processes_mu,lumi_mu,plot_mu["scale_down"],zScaleFac_mu)
		mu_ID=TheStack(processes_mu,lumi_mu,plot_mu["ID"],zScaleFac_mu)
		mu_reso=TheStack(processes_mu,lumi_mu,plot_mu["reso"],zScaleFac_mu)
			
		stackel = TheStack(processes_el,lumi_el,plot_el["default"],zScaleFac_el)
		el_scaleup=TheStack(processes_el,lumi_el,plot_el["scale_up"],zScaleFac_el)
		el_scaledown=TheStack(processes_el,lumi_el,plot_el["scale_down"],zScaleFac_el)
		el_PUup=TheStack(processes_el,lumi_el,plot_el["PU_up"],zScaleFac_el)
		el_PUdown=TheStack(processes_el,lumi_el,plot_el["PU_down"],zScaleFac_el)
	
	if args.znorm:
		muheight = stackmu.theHistogram.FindBin(90)
		print ("Z height of mu: %d +- %d"%(stackmu.theHistogram.GetBinCenter(muheight), stackmu.theHistogram.GetBinWidth(muheight)))
		print ("Z height of mu mc&data: %d, %d"%(stackmu.theHistogram.GetBinContent(muheight), datamu.GetBinContent(muheight)))
		elheight = stackel.theHistogram.FindBin(90)
		print ("Z height of el: %d +- %d"%(stackel.theHistogram.GetBinCenter(elheight), stackel.theHistogram.GetBinWidth(elheight)))
		print ("Z height of el mc&data: %d, %d"%(stackel.theHistogram.GetBinContent(elheight), datael.GetBinContent(elheight)))
		znum = stackmu.theHistogram.GetBinContent(muheight)
		for h in stackmu.theStack.GetHists(): 
			print ("Z height of each hist in mu")
			print (h.GetBinContent(muheight))
			h.Scale(1./znum)
		for h in stackel.theStack.GetHists(): 
			print ("Z height of each hist in ee")
			print (h.GetBinContent(elheight))
			h.Scale(1./znum)
		stackmu.theHistogram.Scale(1./znum)
		stackel.theHistogram.Scale(1./znum)
		
		datamu.Scale(1./znum)
		datael.Scale(1./znum)
	
	if args.ae:
		if args.useall:
			i=0
			for year in range(2016,2019):
				for h in stackmu[i].theStack.GetHists(): inverseAE(h, plot_mu["default"], year)
				for h in stackel[i].theStack.GetHists(): inverseAE(h, plot_el["default"], year)
				inverseAE(stackmu[i].theHistogram, plot_mu["default"], year)
				inverseAE(stackel[i].theHistogram, plot_el["default"], year)
				inverseAE(mu_scaleup[i].theHistogram, plot_mu["scale_up"], year)
				inverseAE(mu_scaledown[i].theHistogram, plot_mu["scale_down"], year)
				inverseAE(mu_ID[i].theHistogram, plot_mu["ID"], year)
				inverseAE(mu_reso[i].theHistogram, plot_mu["reso"], year)
				inverseAE(el_scaleup[i].theHistogram, plot_el["scale_up"], year)
				inverseAE(el_scaledown[i].theHistogram, plot_el["scale_down"], year)
				inverseAE(el_PUup[i].theHistogram, plot_el["PU_up"], year)
				inverseAE(el_PUdown[i].theHistogram, plot_el["PU_down"], year)
				inverseAE(datamu[i], plot_mu["default"], year)
				inverseAE(datael[i], plot_el["default"], year)             
		else:
			if args.use2016: year = 2016
			elif args.use2018: year = 2018
			else: year =2017
			for h in stackmu.theStack.GetHists(): inverseAE(h, plot_mu["default"], year)
			for h in stackel.theStack.GetHists(): inverseAE(h, plot_el["default"], year)
			inverseAE(stackmu.theHistogram, plot_mu["default"], year)
			inverseAE(stackel.theHistogram, plot_el["default"], year)
			inverseAE(mu_scaleup.theHistogram, plot_mu["scale_up"], year)
			inverseAE(mu_scaledown.theHistogram, plot_mu["scale_down"], year)
			inverseAE(mu_ID.theHistogram, plot_mu["ID"], year)
			inverseAE(mu_reso.theHistogram, plot_mu["reso"], year)
			inverseAE(el_scaleup.theHistogram, plot_el["scale_up"], year)
			inverseAE(el_scaledown.theHistogram, plot_el["scale_down"], year)
			inverseAE(el_PUup.theHistogram, plot_el["PU_up"], year)
			inverseAE(el_PUdown.theHistogram, plot_el["PU_down"], year)
			inverseAE(datamu, plot_mu["default"], year)
			inverseAE(datael, plot_el["default"], year)

	if args.useall:
		i=0
		Errs_mu=[]
		Errs_el=[]
		for year in range(2016,2019):
			lis_mu=[[mu_scaleup[i].theHistogram,mu_scaledown[i].theHistogram],mu_ID[i].theHistogram,mu_reso[i].theHistogram]
			lis_el=[[el_scaleup[i].theHistogram,el_scaledown[i].theHistogram],[el_PUup[i].theHistogram,el_PUdown[i].theHistogram]]
			Errs_mu.append(getErrors(stackmu[i].theHistogram,lis_mu))
			Errs_el.append(getErrors(stackel[i].theHistogram,lis_el))
			i+=1
		errmu=Errs_mu[0]+Errs_mu[1]+Errs_mu[2]
		errel=Errs_el[0]+Errs_el[1]+Errs_el[2]         
		stackmu=Addstack(stackmu)
		stackel=Addstack(stackel)
		mu_scaleup=Addstack(mu_scaleup)
		mu_scaledown=Addstack(mu_scaledown)
		mu_ID=Addstack(mu_ID)
		mu_reso=Addstack(mu_reso)
		el_scaleup=Addstack(el_scaleup)
		el_scaledown=Addstack(el_scaledown)
		el_PUup=Addstack(el_PUup)
		el_PUdown=Addstack(el_PUdown)
		datamu=Addhist(datamu)
		datael=Addhist(datael)               
	else:
		lis_mu=[[mu_scaleup.theHistogram,mu_scaledown.theHistogram],mu_ID.theHistogram,mu_reso.theHistogram]
		lis_el=[[el_scaleup.theHistogram,el_scaledown.theHistogram],[el_PUup.theHistogram,el_PUdown.theHistogram]]
		errmu=getErrors(stackmu.theHistogram,lis_mu)
		
		errel=getErrors(stackel.theHistogram,lis_el)
	if args.data:
		yMax = datamu.GetBinContent(datamu.GetMaximumBin())
		if "Mass" in plot_mu["default"].fileName:
			yMin = 0.00001
		else:
			yMin = 0.01
		xMax = datamu.GetXaxis().GetXmax()
		xMin = datael.GetXaxis().GetXmin()
	else:	
		yMax = stackmu.theHistogram.GetBinContent(datahist.GetMaximumBin())
		yMin = 0.01
		xMax = stackmu.theHistogram.GetXaxis().GetXmax()
		xMin = stackmu.theHistogram.GetXaxis().GetXmin()	
	if plot_mu["default"].yMax == None:
		if logScale:
			yMax = yMax*10000
		else:
			yMax = yMax*1.5
	else: yMax = plot_mu["default"].yMax
	
	if "Mass" in plot_mu["default"].fileName:
		yMax = 20000000	
	
	if not plot_mu["default"].yMin == None:
		yMin = plot_mu.yMin
	if not plot_mu["default"].xMin == None:
		xMin = plot_mu["default"].xMin
	if not plot_mu["default"].xMax == None:
		xMax = plot_mu["default"].xMax


	xMin = 200
	xMax = 2000
	yMin = 1e-3
	yMax = 1e4

	if args.ae: 
		yMin = 0.00001 / 40
		yMax = 200000000.0 / 40
		xMin = 200
		xMax = 2000
		yMin *= 10000
		yMax /= 10
	if args.ae:
		# ~ vh = plotPad.DrawFrame(xMin,yMin,xMax,yMax,"; %s ; %s" %("m(l^{+}l^{-}) [GeV]","3 years data"))
		vh = plotPad.DrawFrame(xMin,yMin,xMax,yMax,"; %s ; %s" %("m(l^{+}l^{-}) [GeV]","Events / GeV * 1/(acc. x eff)"))
	else:
		# ~ vh = plotPad.DrawFrame(xMin,yMin,xMax,yMax,"; %s ; %s" %("m(l^{+}l^{-}) [GeV]","Lumi #times d#sigma(pp#rightarrow ll)"))
		vh = plotPad.DrawFrame(xMin,yMin,xMax,yMax,"; %s ; %s" %("m(l^{+}l^{-}) [GeV]","Events / GeV"))
	vh.GetXaxis().SetMoreLogLabels()
	
	drawStack_mu = stackmu
	drawStack_el = stackel

	
	# Draw background from stack
	drawStack_mu.theHistogram.SetLineColor(ROOT.kBlue-3)
	drawStack_el.theHistogram.SetLineColor(ROOT.kRed-3)
	drawStack_mu.theHistogram.Draw("same hist")
	drawStack_el.theHistogram.Draw("same hist")


	# Draw data
	datamu.SetMinimum(0.0001)
	if args.data:
		datamu.SetMarkerColor(ROOT.kViolet)
		datael.SetMarkerColor(ROOT.kOrange)
		datamu.Draw("samep")	
		datael.Draw("samep")

	# Draw legend
	if "Eta" in plot_mu["default"].fileName or "CosTheta" in plot_mu["default"].fileName:
		legendEta.Draw()
	else:
		legend.Draw()

	plotPad.SetLogx(plot_mu["default"].logX)
	if args.useall:
		latex.DrawLatex(0.95, 0.96, "three years data")
	else:	
		latex.DrawLatex(0.95, 0.96, "%.1f fb^{-1} (13 TeV, #mu#mu), %.1f fb^{-1} (13 TeV, ee)"%(lumi_mu*0.001, lumi_el*0.001))
	yLabelPos = 0.85
	cmsExtra = "Private Work"
	if not args.data:
		cmsExtra = "#splitline{Private Work}{Simulation}"
		yLabelPos = 0.82	
	latexCMS.DrawLatex(0.19,0.89,"CMS")
	latexCMSExtra.DrawLatex(0.19,yLabelPos,"%s"%(cmsExtra))
	
	if args.ratio:
		ratioPad.cd()
		ratioPad.SetLogx(plot_mu["default"].logX)

		hhmu = drawStack_mu.theHistogram
		hhel = drawStack_el.theHistogram
		h_mu_scaleup=mu_scaleup.theHistogram
		h_mu_scaledown=mu_scaledown.theHistogram
		h_mu_ID=mu_ID.theHistogram
		h_mu_reso=mu_reso.theHistogram
		h_el_scaleup=el_scaleup.theHistogram
		h_el_scaledown=el_scaledown.theHistogram
		h_el_PUup=el_PUup.theHistogram
		h_el_PUdown=el_PUdown.theHistogram
		ratioGraphs = ROOT.TGraphAsymmErrors(hhmu.GetSize()-2)
		chann = True if "BB" in plot_mu["default"].fileName else False
		#print(errel)
		#print(errmu)
		for i in range(1, hhmu.GetSize()-1):
			xval = hhmu.GetBinCenter(i)
			xerr = hhmu.GetBinWidth(i)/2
			if hhel.GetBinContent(i) == 0: continue
			if hhmu.GetBinContent(i) == 0: continue
			#print(math.sqrt(errmu[i-1])/hhmu.GetBinContent(i))
			yval = hhmu.GetBinContent(i)*1.0/hhel.GetBinContent(i)
			yerr = yval*math.sqrt((errel[i-1]**0.5/hhel.GetBinContent(i))**2+(errmu[i-1]**0.5/hhmu.GetBinContent(i))**2+(hhmu.GetBinError(i)/hhmu.GetBinContent(i))**2+(hhel.GetBinError(i)/hhel.GetBinContent(i))**2)
			ratioGraphs.SetPoint(i, xval, yval)
			ratioGraphs.SetPointError(i, xerr, xerr, yerr, yerr)
			print ("M = %f, r+-e = %f +- %f"%(xval, yval, yerr/yval))
		ratioData = ROOT.TGraphAsymmErrors(datamu.GetSize()-2)
		for i in range(1, datamu.GetSize()-1):
			xval = datamu.GetBinCenter(i)
			xerr = datamu.GetBinWidth(i)/2
			if datael.GetBinContent(i) == 0: continue
			if datamu.GetBinContent(i) == 0: continue
			#print(datael.GetBinContent(i))
			#print(datael.GetBinError(i))
			yval = datamu.GetBinContent(i)*1.0/datael.GetBinContent(i)
			yerr = yval*math.sqrt((datamu.GetBinError(i)/datamu.GetBinContent(i))**2+(datael.GetBinError(i)/datael.GetBinContent(i))**2)
			ratioData.SetPoint(i, xval, yval)
			ratioData.SetPointError(i, xerr, xerr, yerr, yerr)
			print ("M = %f, r+-e = %f +- %f"%(xval, yval, yerr/yval))
		nBinsX = 20
		nBinsY = 10
		if args.ae: 
			hAxis = ROOT.TH2F("hAxis", "", nBinsX, xMin, xMax, nBinsY, 0.5, 2.5)
		else:	
			hAxis = ROOT.TH2F("hAxis", "", nBinsX, xMin, xMax, nBinsY, 0.5, 5)
		hAxis.Draw("AXIS")

		hAxis.GetYaxis().SetNdivisions(408)
		hAxis.SetTitleOffset(0.4, "Y")
		hAxis.SetTitleSize(0.09, "Y")
		hAxis.SetTitleSize(0.06, "X")
		hAxis.SetYTitle("R_{#mu#mu/ee}")
		hAxis.SetXTitle("m(l^{+}l^{-}) [GeV]")
		hAxis.GetXaxis().SetLabelSize(0.048)
		hAxis.GetYaxis().SetLabelSize(0.048)
		#hAxis.GetXaxis().SetTicks("+")
		#hAxis.SetTitleSize(0.15, "Y")
		hAxis.GetXaxis().SetMoreLogLabels()	
		oneLine = ROOT.TLine(xMin, 1.0, xMax, 1.0)
		oneLine.SetLineStyle(2)
		oneLine.Draw()
	
		ratioGraphs.SetFillColor(ROOT.kBlue-3)
		ratioGraphs.SetMarkerColor(ROOT.kBlue-3)
		ratioGraphs.GetXaxis().SetLabelSize(0.0)
		ratioGraphs.SetFillStyle(3002)	
		ratioGraphs.Draw("SAME p")
		ratioData.SetMarkerColor(ROOT.kViolet)
		ratioData.Draw("same p")
		
		rlegend = TLegend(0.2, 0.65, 0.5, 0.925)
		rlegend.SetFillStyle(0)
		rlegend.SetBorderSize(1)
		rlegend.SetTextFont(42)
		rlegend.AddEntry(ratioGraphs, "e/mu in MC inclusive", "pe")
		rlegend.AddEntry(ratioData, "e/mu in data", "pe")
		rlegend.Draw("same")
		
		ratioPad.Update()
					

	ROOT.gPad.RedrawAxis()
	plotPad.RedrawAxis()
	if args.ratio:
		ratioPad.RedrawAxis()
	if not os.path.exists("lepFlavor"):
		os.makedirs("lepFlavor")	

	if args.use2016: year = "2016"
	elif args.useall: year = "2016_to_2018"
	elif args.use2018: year = "2018"
	else: year = "2017"

	if args.ae: year += "_inverseAE"
	if args.znorm: year += "_znorm"
	
	hCanvas.Print("lepFlavor/%s_%s_datamc.pdf"%(plot_mu["default"].fileName, year))


					
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
	parser.add_argument("-b", "--backgrounds", dest="backgrounds", action="append", default=[],
						  help="backgrounds to plot.")
	parser.add_argument("--ae", action="store_true", dest="ae", default=False,help="times inverse Acceptance x Efficiency")
	parser.add_argument("--znorm", action="store_true", dest="znorm", default=False, help="normalize to z peak")
	parser.add_argument("--all", action="store_true", dest="useall", default=False, help="add the data from 2016 to 2018")
	args = parser.parse_args()
	if len(args.backgrounds) == 0:
		args.backgrounds = ["Wjets","Other","DrellYan"]

		
	muplots_bb={"default":"massPlotBB","scale_up":"massPlotBBScaleUpNoLog","scale_down":"massPlotBBScaleDownNoLog","reso":"massPlotBBSmearNoLog","ID":"massPlotBBMuonIDNoLog"}
	muplots_be={"default":"massPlotBE","scale_up":"massPlotBEScaleUpNoLog","scale_down":"massPlotBEScaleDownNoLog","reso":"massPlotBESmearNoLog","ID":"massPlotBEMuonIDNoLog"}
	eleplots_bb={"default":"massPlotEleBB","scale_up":"massPlotEleBBScaleUpNoLog","scale_down":"massPlotEleBBScaleDownNoLog","PU_up":"massPlotEleBBPUScaleUpNoLog","PU_down":"massPlotEleBBPUScaleDownNoLog"}
	eleplots_be={"default":"massPlotEleBE","scale_up":"massPlotEleBEScaleUpNoLog","scale_down":"massPlotEleBEScaleDownNoLog","PU_up":"massPlotEleBEPUScaleUpNoLog","PU_down":"massPlotEleBEPUScaleDownNoLog"}
	plot_mu_bb={}
	plot_el_bb={}
	plot_mu_be={}
	plot_el_be={}
	for key in muplots_bb.keys():
		plot_mu_bb[key] = getPlot(muplots_bb[key])
		plot_mu_bb[key].logX=True
	for key in eleplots_bb.keys():
		plot_el_bb[key] = getPlot(eleplots_bb[key])
		plot_el_bb[key].logX=True
	for key in muplots_be.keys():
		plot_mu_be[key] = getPlot(muplots_be[key])
		plot_mu_be[key].logX=True
	for key in eleplots_be.keys():
		plot_el_be[key] = getPlot(eleplots_be[key])
		plot_el_be[key].logX=True

	plotDataMC(args,plot_mu_bb,plot_el_bb)
	plotDataMC(args,plot_mu_be,plot_el_be)
