from ROOT import * 
from numpy import array as ar
from array import array
#from setTDRStyle import setTDRStyle
from copy import deepcopy
import pickle
from helpers import *
from defs import getPlot, Backgrounds, Signals, Data


tableTemplate = '''
\\begin{table}
\\begin{center}
\\begin{tabular}{c|ccccc}
 Signal Model & \multicolumn{5}{c}{BB category} \\\\ \\hline
 & 400-500 & 500-700 & 700-1100 & 1100-1900 & 1900-3500 \\\\ \\hline
 %s 
 %s
 %s
 %s
 %s
 %s 
 & \multicolumn{5}{c}{BE category} \\\\ \\hline
 & 400-500 & 500-700 & 700-1100 & 1100-1900 & 1900-3500 \\\\ \\hline
 %s 
 %s
 %s
 %s
 %s
 %s 
\end{tabular}
\end{center}
\end{table}
'''
lineTemplate = " %s & %.2f \\pm %.2f  & %.2f \\pm %.2f & %.2f \\pm %.2f & %.2f \\pm %.2f & %.2f \\pm %.2f  \\\\"


def main():
	gROOT.SetBatch(True)

	
	histos = ["BB","BE"]
	labels = ["dimuon_BB","dimuon_BE"]
	
	lambdas = [10,16,22,28,34]
	models = ["ConLL","ConLR","ConRR"]
	bins = [4,7]

	massPlot = getPlot("massPlotForLimit")
	massPlotSmeared = getPlot("massPlotSmeared")
	massPlotUp = getPlot("massPlotUp")
	massPlotDown = getPlot("massPlotDown")
	massPlotWeighted = getPlot("massPlotWeighted")
	xMax = 5000
	massBins = [400,500,700,1100,1900,3500,xMax]
	massCuts = [400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,2600,2700,2800,2900,3000,3100,3200,3300,3400,3500]
	
	signalYields = {}
	signalYieldsSingleBin = {}
	
	names = ["default","resolution","scale","ID"]
	
	
	
	massPlots = [massPlot,massPlotSmeared,massPlotUp,massPlotWeighted]
	
	for number, plot in enumerate(massPlots):
		suffix = names[number]
		for i, histo in enumerate(histos):
			label = labels[i]

			for model in models:
				for l in lambdas:
					name = "CITo2Mu_Lam%dTeV%s"%(l,model)
					signal = Process(getattr(Signals,name))
					signalYields["%s_%s"%(name,label)] = {}
					signalYieldsSingleBin["%s_%s"%(name,label)] = {}
					nameDY = "CITo2Mu_Lam100kTeV%s"%(model)
					signalDY = Process(getattr(Signals,nameDY))
					
					sigHist = deepcopy(signal.loadHistogramProjected(plot, bins[i]))
					if "_BE" in label:
						sigHist.Add(deepcopy(signal.loadHistogramProjected(plot, 10)))					
						
					sigHistDY = deepcopy(signalDY.loadHistogramProjected(plot, bins[i]))
					if "_BE" in label:
						sigHistDY.Add(deepcopy(signalDY.loadHistogramProjected(plot, 10)))
		
					xMin = 300
					if "Des" in name:
						xMin = 300
		
					func = TF1("f1_%s"%name,"[5]*TMath::Exp([0]+[1]*x+[2]*x^2+[3]*x^3)*x^[4]",xMin,5000)
					func.SetParameter(0,30)
					func.SetParameter(1,-0.001)
					func.SetParameter(2,2e-07)
					func.SetParameter(3,-2e-11)
					func.SetParameter(4,-5)
					func.SetParameter(5,100)

					funcDY = TF1("f2_%s"%name,"[5]*TMath::Exp([0]+[1]*x+[2]*x^2+[3]*x^3)*x^[4]",xMin,5000)
					funcDY.SetParameter(0,30)
					funcDY.SetParameter(1,-0.001)
					funcDY.SetParameter(2,2e-07)
					funcDY.SetParameter(3,-2e-11)
					funcDY.SetParameter(4,-5)
					funcDY.SetParameter(5,100)
					
					
					
		
					hCanvas = TCanvas("hCanvas", "Distribution", 800,800)
					
					plotPad = ROOT.TPad("plotPad","plotPad",0,0,1,1)
#					setTDRStyle()		
#					plotPad.UseCurrentStyle()
					plotPad.Draw()	
					plotPad.cd()
					sigHist.Rebin(20)
					sigHist.Fit("f1_%s"%name,"MR")				
		
					for index, mass in enumerate(massBins):
						if index < len(massBins)-1:
							valueSig = func.Integral(mass,massBins[index+1])/20
							errorSig = func.IntegralError(mass,massBins[index+1])/20
							signalYields["%s_%s"%(name,label)][str(index)] = [valueSig,errorSig]
					for index, mass in enumerate(massCuts):
						valueSig = func.Integral(mass,5000)/20
						errorSig = func.IntegralError(mass,5000)/20
						signalYieldsSingleBin["%s_%s"%(name,label)][str(mass)] = [valueSig,errorSig]
					sigHistDY.Rebin(20)
					sigHistDY.Fit("f2_%s"%name,"MR")		
											
					for index, mass in enumerate(massBins):
						if index < len(massBins)-1:
							valueDY = funcDY.Integral(mass,massBins[index+1])/20
							errorDY = funcDY.IntegralError(mass,massBins[index+1])/20
							diff = signalYields["%s_%s"%(name,label)][str(index)][0] - valueDY
							errDiff = diff*((signalYields["%s_%s"%(name,label)][str(index)][1]/signalYields["%s_%s"%(name,label)][str(index)][0])**2 + (errorDY/valueDY)**2)**0.5
							signalYields["%s_%s"%(name,label)][str(index)] = [diff,errDiff]
					for index, mass in enumerate(massCuts):
						valueDY = funcDY.Integral(mass,5000)/20
						errorDY = funcDY.IntegralError(mass,5000)/20
						diff = signalYieldsSingleBin["%s_%s"%(name,label)][str(mass)][0] - valueDY
						errDiff = diff*((signalYieldsSingleBin["%s_%s"%(name,label)][str(mass)][1]/signalYieldsSingleBin["%s_%s"%(name,label)][str(mass)][0])**2 + (errorDY/valueDY)**2)**0.5
						signalYieldsSingleBin["%s_%s"%(name,label)][str(mass)] = [diff,errDiff]
							

					plotPad.DrawFrame(xMin,0.001,3500,3500,"; m [GeV]; Events")	
					plotPad.SetLogy()			
					
					
					sigHist.Draw("samehist")
					sigHistDY.Draw("samehist")
					sigHistDY.SetLineColor(kBlack)
					func.SetLineColor(sigHist.GetLineColor())
					funcDY.SetLineColor(sigHistDY.GetLineColor())
					func.Draw("same")
					funcDY.Draw("same")
					
					
					legend = TLegend(0.375, 0.7, 0.925, 0.925)
					legend.SetFillStyle(0)
					legend.SetBorderSize(0)
					legend.SetTextFont(42)		
					legend.AddEntry(sigHist,"%s"%name,"l")	
					legend.AddEntry(sigHistDY,"%s"%nameDY,"l")	
		
					legend.Draw()
						
					
					hCanvas.Print("plots/fit_%s_%s_%s.pdf"%(name,label,suffix))



		outFilePkl = open("signalYields_%s.pkl"%suffix,"w")
		pickle.dump(signalYields, outFilePkl)
		outFilePkl.close()		
		outFilePkl = open("signalYieldsSingleBin_%s.pkl"%suffix,"w")
		pickle.dump(signalYieldsSingleBin, outFilePkl)
		outFilePkl.close()		
	
			
							
main()
