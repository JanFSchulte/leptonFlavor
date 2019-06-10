from ROOT import *

from setTDRStyle import setTDRStyle

dyFile = TFile("drellyan.root","OPEN")	
dyHist = dyFile.Get("h1_dimuo_m")
dyHistPDF = dyFile.Get("h1_dimuo_m_pdf_total")


otherFile = TFile("other.root","OPEN")	
otherHist = otherFile.Get("h1_dimuo_m")
otherHistPDF = otherFile.Get("h1_dimuo_m_pdf_total")

from setTDRStyle import setTDRStyle

xMin = 2000
xMax = 60000

defaultDY =  dyHist.Integral(dyHist.FindBin(xMin),dyHist.FindBin(xMax-0.01))
pdfDY     =  dyHistPDF.Integral(dyHistPDF.FindBin(xMin),dyHistPDF.FindBin(xMax-0.01))

defaultOther =  otherHist.Integral(otherHist.FindBin(xMin),otherHist.FindBin(xMax-0.01))
pdfOther     =  otherHistPDF.Integral(otherHistPDF.FindBin(xMin),otherHistPDF.FindBin(xMax-0.01))

print "DY PDF: ", abs(1 - (pdfDY+pdfOther)/(defaultDY+defaultOther))

xMin = 4000
xMax = 60000

defaultDY =  dyHist.Integral(dyHist.FindBin(xMin),dyHist.FindBin(xMax-0.01))
pdfDY     =  dyHistPDF.Integral(dyHistPDF.FindBin(xMin),dyHistPDF.FindBin(xMax-0.01))

defaultOther =  otherHist.Integral(otherHist.FindBin(xMin),otherHist.FindBin(xMax-0.01))
pdfOther     =  otherHistPDF.Integral(otherHistPDF.FindBin(xMin),otherHistPDF.FindBin(xMax-0.01))

print "DY PDF: ", abs(1 - (pdfDY+pdfOther)/(defaultDY+defaultOther))

massBins = [400,500,700,1100,1900,3500,5000]
pdfUncertDY = []
pdfUncertOther = []
for i in range(0,len(massBins)-1):
	xMin = massBins[i]
	xMax = massBins[i+1]
	print xMin, xMax
	defaultDY =  dyHist.Integral(dyHist.FindBin(xMin),dyHist.FindBin(xMax-0.01))
	pdfDY     =  dyHistPDF.Integral(dyHistPDF.FindBin(xMin),dyHistPDF.FindBin(xMax-0.01))
	print defaultDY, pdfDY
	defaultOther =  otherHist.Integral(otherHist.FindBin(xMin),otherHist.FindBin(xMax-0.01))
	pdfOther     =  otherHistPDF.Integral(otherHistPDF.FindBin(xMin),otherHistPDF.FindBin(xMax-0.01))

	pdfUncertDY.append(1+abs(1 - pdfDY/defaultDY))
	pdfUncertOther.append(1+abs(1 - pdfOther/defaultOther))
	
print 	pdfUncertDY
print  pdfUncertOther


hCanvas = TCanvas("hCanvas", "Distribution", 800,400)

plotPad = TPad("plotPad","plotPad",0,0,1,1)
setTDRStyle()		
plotPad.UseCurrentStyle()
plotPad.Draw()	
plotPad.cd()


dyHist.Rebin(500)
dyHistPDF.Rebin(500)
otherHist.Rebin(500)
otherHistPDF.Rebin(500)
#~ 
#~ 
dyHistPDF.Divide(dyHist)
otherHistPDF.Divide(otherHist)


plotPad.DrawFrame(0,0.95,5000,1.5,"; dilepton mass [GeV]; PDF uncertainty")

dyHistPDF.SetLineColor(kGreen+2)
otherHistPDF.SetLineColor(kRed)
dyHist.GetXaxis().SetRangeUser(0,2000)
dyHistPDF.Draw("samehist")
#~ dyHistPDF.Draw("samepe")
otherHistPDF.Draw("samehist")

legend = TLegend(0.2, 0.7, 0.5, 0.925)
legend.SetFillStyle(0)
legend.SetBorderSize(0)
legend.SetTextFont(42)		
legend.AddEntry(dyHistPDF,"Drell-Yan","l")		
legend.AddEntry(otherHistPDF,"Other Backgrounds","l")		
legend.Draw()



	
hCanvas.Print("pdfUncert.pdf")
