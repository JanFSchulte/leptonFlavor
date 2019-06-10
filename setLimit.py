import ROOT
from ROOT import TFile, TTree, TCanvas, TGraph, TMultiGraph, TGraphErrors, TLegend
import CMS_lumi, tdrstyle
import subprocess # to execute shell command
ROOT.gROOT.SetBatch(ROOT.kTRUE)
 
# CMS style
CMS_lumi.cmsText = "CMS"
CMS_lumi.extraText = "Preliminary"
CMS_lumi.cmsTextSize = 0.65
CMS_lumi.outOfFrame = True
tdrstyle.setTDRStyle()
 
 
# CREATE datacards
def createDataCardsThetaB(masses, ysigs, ybkgs, yerrs):
 
    datacard_lines = """# automatic generated counting experiment
                       imax 1  number of channels
                       jmax 1  number of backgrounds
                       kmax 1  number of nuisance parameters (sources of systematical uncertainties)
                       ------------
                       bin b1
                       observation 0
                       ------------
                       bin              b1     b1    
                       process         sig   bkg  
                       process          0     1    
                       rate            %f  %f  
                       ------------
                       lumi    lnN    %f  %f   uncert
                       """
 
    # make datacards for differents values of theta_B
    for i in range(len(masses)):
        datacard = open("datacard_"+str(int(masses[i]))+".txt", 'w')
        datacard.write(datacard_lines%(ysigs[i], ybkgs[i], 1+yerrs[i], 1+yerrs[i]))
        datacard.close()
        #print ">>>   datacard_"+label+".txt created."
 
 
# EXECUTE datacards
def executeDataCards(masses):
 
    for mass in masses:
        file_name = "datacard_"+str(int(mass))+".txt"
        combine_command = "combine -M Asymptotic -m %s %s" % (int(mass),file_name)
        print ""
        print ">>> " + combine_command
        p = subprocess.Popen(combine_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in p.stdout.readlines():
            print line.rstrip("\n")
        #print ">>>   higgsCombine"+label+".Asymptotic.mH125.root created"
        retval = p.wait()
 
 
# GET limits from root file
def getLimits(file_name):
 
    file = TFile(file_name)
    tree = file.Get("limit")
 
    limits = [ ]
    for quantile in tree:
        limits.append(tree.limit)
        print ">>>   %.2f" % limits[-1]
 
    return limits[:6]
 
 
# PLOT upper limits
def plotUpperLimits(masses):
    # see CMS plot guidelines: https://ghm.web.cern.ch/ghm/plots/
 
    N = len(masses)
    yellow = TGraph(2*N)    # yellow band
    green = TGraph(2*N)     # green band
    median = TGraph(N)      # median line
    obs = TGraph(N)
 
    up2s = [ ]
    for i in range(N):
        file_name = "higgsCombineTest.Asymptotic.mH%d.root"%int(masses[i])
        limit = getLimits(file_name)
        up2s.append(limit[4])
        yellow.SetPoint(    i,    masses[i], limit[4] ) # + 2 sigma
        green.SetPoint(     i,    masses[i], limit[3] ) # + 1 sigma
        median.SetPoint(    i,    masses[i], limit[2] ) # median
        green.SetPoint(  2*N-1-i, masses[i], limit[1] ) # - 1 sigma
        yellow.SetPoint( 2*N-1-i, masses[i], limit[0] ) # - 2 sigma
        obs.SetPoint(i, masses[i], limit[5]) 

    W = 800
    H  = 600
    T = 0.08*H
    B = 0.12*H
    L = 0.12*W
    R = 0.04*W
    c = TCanvas("c","c",100,100,W,H)
    c.SetFillColor(0)
    c.SetBorderMode(0)
    c.SetFrameFillStyle(0)
    c.SetFrameBorderMode(0)
    c.SetLeftMargin( L/W )
    c.SetRightMargin( R/W )
    c.SetTopMargin( T/H )
    c.SetBottomMargin( B/H )
    c.SetTickx(0)
    c.SetTicky(0)
    c.SetGrid()
    c.cd()
    frame = c.DrawFrame(1.4,0.001, 4.1, 10)
    frame.GetYaxis().CenterTitle()
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetTitleOffset(0.9)
    frame.GetXaxis().SetNdivisions(508)
    frame.GetYaxis().CenterTitle(True)
    frame.GetYaxis().SetTitle("95% upper limit on #sigma / #sigma_{SM}")
#    frame.GetYaxis().SetTitle("95% upper limit on #sigma #times BR / (#sigma #times BR)_{SM}")
    frame.GetXaxis().SetTitle("M_{ll} [GeV]")
    frame.SetMinimum(0)
    frame.SetMaximum(max(up2s)*1.5)
    frame.GetXaxis().SetLimits(min(masses),max(masses))
 
    yellow.SetFillColor(ROOT.kOrange)
    yellow.SetLineColor(ROOT.kOrange)
    yellow.SetFillStyle(1001)
    yellow.Draw('F')
 
    green.SetFillColor(ROOT.kGreen+1)
    green.SetLineColor(ROOT.kGreen+1)
    green.SetFillStyle(1001)
    green.Draw('Fsame')
 
    median.SetLineColor(1)
    median.SetLineWidth(2)
    median.SetLineStyle(2)
    median.Draw('Lsame')
    
    obs.SetLineColor(2)
    obs.SetLineWidth(3)
    #obs.SetLineStyle(2)
    obs.Draw("Lsame")
 
    CMS_lumi.CMS_lumi(c,13,11)
    ROOT.gPad.SetTicks(1,1)
    frame.Draw('sameaxis')
 
    x1 = 0.15
    x2 = x1 + 0.24
    y2 = 0.85
    y1 = 0.69
    legend = TLegend(x1,y1,x2,y2)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.041)
    legend.SetTextFont(42)
    legend.AddEntry(obs, "Observed Limit (BB)", "L")
    legend.AddEntry(median, "Asymptotic CL_{s} expected",'L')
    legend.AddEntry(green, "#pm 1 std. deviation",'f')
#    legend.AddEntry(green, "Asymptotic CL_{s} #pm 1 std. deviation",'f')
    legend.AddEntry(yellow,"#pm 2 std. deviation",'f')
#    legend.AddEntry(green, "Asymptotic CL_{s} #pm 2 std. deviation",'f')
    legend.Draw()
 
    print " "
    c.SaveAs("LFU_Limit_BB.pdf")
    c.Close()
 
 
# RANGE of floats
def frange(start, stop, step):
    i = start
    while i <= stop:
        yield i
        i += step
 
 
# MAIN
def main():
 
    masses = [1081, 1170, 1267, 1371, 1485, 1607, 1740, 1883]
    ysigs = [1.017, 1.059, 1.105, 1.090, 1.057, 1.242, 1.597, 1.009]
    ybkgs = [1.089, 1.297, 1.380, 1.259, 1.024, 0.451, 1.344, 0.999]
    yerrs = [0.083, 0.083, 0.083, 0.083, 0.083, 0.083, 0.083, 0.083]
    #ysigs = [0.918, 0.905, 1.187, 0.861, 1.162, 1.111, 1.027, 1.212]
    #ybkgs = [1.019, 1.045, 0.851, 0.756, 1.194, 1.496, 1.053, 0.415]
    #yerrs = [0.094, 0.094, 0.094, 0.094, 0.094, 0.094, 0.094, 0.094]
 
    createDataCardsThetaB(masses,ysigs,ybkgs,yerrs)
    executeDataCards(masses)
    plotUpperLimits(masses)
 
 
 
if __name__ == '__main__':
    main()
