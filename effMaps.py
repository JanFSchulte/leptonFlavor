from ROOT import * 
from numpy import array as ar
from array import array
from setTDRStyle import setTDRStyle
from copy import deepcopy

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
lineTemplate = " %s & %.3f  & %.3f & %.3f & %.3f & %.3f  \\\\"

def main():
	### for data

	


	genPlot = getPlot("etaPtMapGen")
	recoPlot = getPlot("etaPtMapReco")


	drellyanGen = Process(getattr(Backgrounds,"DrellYan"))
	drellyanReco = Process(getattr(Backgrounds,"DrellYan"))

	genPlot = deepcopy(drellyanGen.loadHistogram(genPlot))
	recoPlot = deepcopy(drellyanReco.loadHistogram(recoPlot))

	

	eff = recoPlot.Clone("bla")
	eff.Divide(genPlot)

	print genPlot.GetBinContent(3,5)
	print recoPlot.GetBinContent(3,5)

	hCanvas = TCanvas("hCanvas", "Distribution", 800,800)

	plotPad = ROOT.TPad("plotPad","plotPad",0,0,1,1)

	style = setTDRStyle()		
	gStyle.SetOptStat(0)
	gStyle.SetPadRightMargin(0.2)	
	gStyle.SetPadBottomMargin(0.1)	
	plotPad.UseCurrentStyle()
	gStyle.SetOptStat(0)
	gStyle.SetPadRightMargin(0.2)
	plotPad.Draw()	
	plotPad.cd()
	
	#~ plotPad.DrawFrame(-2.4,0,2.4,2000,"; muon #eta ; muon p_{T} GeV")
	eff.GetZaxis().SetRangeUser(0,1)
	eff.GetYaxis().SetRangeUser(53,2000)
	eff.GetZaxis().SetTitle("Trigger+Reconstruction+ID efficiency")
	eff.GetZaxis().SetTitleOffset(1.3)
	eff.GetYaxis().SetTitle("simulated muon p_{T} GeV [GeV]")
	eff.GetYaxis().SetTitleOffset(2)
	eff.GetXaxis().SetTitle("simulated muon #eta")
	eff.Draw("colz")
	plotPad.SetLogy()

	latex = TLatex()
	latex.SetTextFont(42)
	latex.SetTextAlign(31)
	latex.SetTextSize(0.04)
	latex.SetNDC(True)
	latexCMS = TLatex()
	latexCMS.SetTextFont(61)
	latexCMS.SetTextSize(0.055)
	latexCMS.SetNDC(True)
	latexCMSExtra = TLatex()
	latexCMSExtra.SetTextFont(52)
	latexCMSExtra.SetTextSize(0.03)
	latexCMSExtra.SetNDC(True) 
		
	latex.DrawLatex(0.85, 0.96, "(13 TeV)")
	
	cmsExtra = "Supplementary"
	latexCMS.DrawLatex(0.15,0.96,"CMS")
	if "Simulation" in cmsExtra:
		yLabelPos = 0.81	
	else:
		yLabelPos = 0.96	

	latexCMSExtra.DrawLatex(0.27,yLabelPos,"%s"%(cmsExtra))	

	hCanvas.Print("effMapMuons.pdf")
	
	

main()
