from ROOT import *
from setTDRStyle import setTDRStyle
import pickle
pkl = open("signalYields_default.pkl", "r")
signalYields_default = pickle.load(pkl)
pkl = open("signalYields_scale.pkl", "r")
signalYields_scale = pickle.load(pkl)

pkl = open("signalYields_resolution.pkl", "r")
signalYields_resolution = pickle.load(pkl)

pkl = open("signalYields_ID.pkl", "r")
signalYields_ID = pickle.load(pkl)

histos = ["BB","BE"]
labels = ["dimuon_Moriond2017_BB","dimuon_Moriond2017_BE"]
lambdas = [10,16,22,28,34]
models = ["ConLL","ConLR","ConRR","DesLL","DesLR","DesRR"]

scaleHist1 = TH1F("scaleHist1","scaleHist1",25,0.9,1.1)
resolutionHist1 = TH1F("smearHist1","smearHist1",25,0.9,1.1)
IDHist1 = TH1F("IDHist1","IDHist1",25,0.9,1.1)
scaleHist2 = TH1F("scaleHist2","scaleHist2",25,0.9,1.1)
resolutionHist2 = TH1F("smearHist2","smearHist2",25,0.9,1.1)
IDHist2 = TH1F("IDHist2","IDHist2",25,0.9,1.1)
scaleHist3 = TH1F("scaleHist3","scaleHist3",25,0.9,1.1)
resolutionHist3 = TH1F("smearHist3","smearHist3",25,0.9,1.1)
IDHist3 = TH1F("IDHist3","IDHist3",25,0.9,1.1)
scaleHist4 = TH1F("scaleHist4","scaleHist4",25,0.9,1.1)
resolutionHist4 = TH1F("smearHist4","smearHist4",25,0.9,1.1)
IDHist4 = TH1F("IDHist4","IDHist4",25,0.9,1.1)
scaleHist5 = TH1F("scaleHist5","scaleHist5",25,0.9,1.1)
resolutionHist5 = TH1F("smearHist5","smearHist5",25,0.9,1.1)
IDHist5 = TH1F("IDHist5","IDHist5",25,0.9,1.1)

for i, histo in enumerate(histos):
	label = labels[i]

	for model in models:
		for l in lambdas:
			name = "CITo2Mu_Lam%dTeV%s_%s"%(l,model,label)
			
			print signalYields_scale[name][str(0)][0]/signalYields_default[name][str(0)][0]
			
			scaleHist1.Fill(signalYields_scale[name][str(0)][0]/signalYields_default[name][str(0)][0])
			scaleHist2.Fill(signalYields_scale[name][str(1)][0]/signalYields_default[name][str(1)][0])
			scaleHist3.Fill(signalYields_scale[name][str(2)][0]/signalYields_default[name][str(2)][0])
			scaleHist4.Fill(signalYields_scale[name][str(3)][0]/signalYields_default[name][str(3)][0])
			scaleHist5.Fill(signalYields_scale[name][str(4)][0]/signalYields_default[name][str(4)][0])
			
			resolutionHist1.Fill(signalYields_resolution[name][str(0)][0]/signalYields_default[name][str(0)][0])
			resolutionHist2.Fill(signalYields_resolution[name][str(1)][0]/signalYields_default[name][str(1)][0])
			resolutionHist3.Fill(signalYields_resolution[name][str(2)][0]/signalYields_default[name][str(2)][0])
			resolutionHist4.Fill(signalYields_resolution[name][str(3)][0]/signalYields_default[name][str(3)][0])
			resolutionHist5.Fill(signalYields_resolution[name][str(4)][0]/signalYields_default[name][str(4)][0])

			IDHist1.Fill(signalYields_ID[name][str(0)][0]/signalYields_default[name][str(0)][0])
			IDHist2.Fill(signalYields_ID[name][str(1)][0]/signalYields_default[name][str(1)][0])
			IDHist3.Fill(signalYields_ID[name][str(2)][0]/signalYields_default[name][str(2)][0])
			IDHist4.Fill(signalYields_ID[name][str(3)][0]/signalYields_default[name][str(3)][0])
			IDHist5.Fill(signalYields_ID[name][str(4)][0]/signalYields_default[name][str(4)][0])
			
scaleHist1.SetLineColor(kRed)			
scaleHist2.SetLineColor(kBlue)			
scaleHist3.SetLineColor(kBlack)			
scaleHist4.SetLineColor(kOrange+1)			
scaleHist5.SetLineColor(kGreen+3)			

resolutionHist1.SetLineColor(kRed)			
resolutionHist2.SetLineColor(kBlue)			
resolutionHist3.SetLineColor(kBlack)			
resolutionHist4.SetLineColor(kOrange+1)			
resolutionHist5.SetLineColor(kGreen+3)			

IDHist1.SetLineColor(kRed)			
IDHist2.SetLineColor(kBlue)			
IDHist3.SetLineColor(kBlack)			
IDHist4.SetLineColor(kOrange+1)			
IDHist5.SetLineColor(kGreen+3)			

hCanvas = TCanvas("hCanvas", "Distribution", 800,800)	

plotPad = TPad("plotPad","plotPad",0,0,1,1)
#~ ratioPad = ROOT.TPad("ratioPad","ratioPad",0,0.,1,0.3)
setTDRStyle()		
plotPad.UseCurrentStyle()
#~ ratioPad.UseCurrentStyle()
plotPad.Draw()	
#~ ratioPad.Draw()	
plotPad.cd()

plotPad.DrawFrame(0.9,0,1.1,20,"; impact of scale systematic; Signal models")

scaleHist1.Draw("samehist")
scaleHist2.Draw("samehist")
scaleHist3.Draw("samehist")
scaleHist4.Draw("samehist")
scaleHist5.Draw("samehist")


legend = TLegend(0.5, 0.7, 0.925, 0.925)
legend.SetFillStyle(0)
legend.SetBorderSize(0)
legend.SetTextFont(42)		
legend.AddEntry(scaleHist1,"400-500 GeV","l")	
legend.AddEntry(scaleHist2,"500-700 GeV","l")	
legend.AddEntry(scaleHist3,"700-1100 GeV","l")	
legend.AddEntry(scaleHist4,"1100-1900 GeV","l")	
legend.AddEntry(scaleHist5,"1900-3500 GeV","l")	
legend.Draw()

hCanvas.Print("scaleUncert.pdf")

plotPad.DrawFrame(0.9,0,1.1,30,"; impact of resolution systematic; Signal models")

resolutionHist1.Draw("samehist")
resolutionHist2.Draw("samehist")
resolutionHist3.Draw("samehist")
resolutionHist4.Draw("samehist")
resolutionHist5.Draw("samehist")


legend = TLegend(0.5, 0.7, 0.925, 0.925)
legend.SetFillStyle(0)
legend.SetBorderSize(0)
legend.SetTextFont(42)		
legend.AddEntry(resolutionHist1,"400-500 GeV","l")	
legend.AddEntry(resolutionHist2,"500-700 GeV","l")	
legend.AddEntry(resolutionHist3,"700-1100 GeV","l")	
legend.AddEntry(resolutionHist4,"1100-1900 GeV","l")	
legend.AddEntry(resolutionHist5,"1900-3500 GeV","l")	
legend.Draw()

hCanvas.Print("resolutionUncert.pdf")
plotPad.DrawFrame(0.9,0,1.1,30,"; impact of ID systematic; Signal models")

IDHist1.Draw("samehist")
IDHist2.Draw("samehist")
IDHist3.Draw("samehist")
IDHist4.Draw("samehist")
IDHist5.Draw("samehist")


legend = TLegend(0.5, 0.7, 0.925, 0.925)
legend.SetFillStyle(0)
legend.SetBorderSize(0)
legend.SetTextFont(42)		
legend.AddEntry(IDHist1,"400-500 GeV","l")	
legend.AddEntry(IDHist2,"500-700 GeV","l")	
legend.AddEntry(IDHist3,"700-1100 GeV","l")	
legend.AddEntry(IDHist4,"1100-1900 GeV","l")	
legend.AddEntry(IDHist5,"1900-3500 GeV","l")	
legend.Draw()

hCanvas.Print("IDUncert.pdf")
