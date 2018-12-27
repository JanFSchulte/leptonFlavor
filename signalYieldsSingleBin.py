from ROOT import gROOT, TFile
from numpy import array as ar
from array import array
from copy import deepcopy
import pickle



def main():
	gROOT.SetBatch(True)

	
	histos = ["BB","BE"]
	labels = ["dimuon_2017","dielectron_2017"]
	#~ channels = ["cito2mu","cito2e"]
	suffixesMu = ["nominal","scaledown","smeared","muonid"]
	suffixesEle = ["nominal","scaledown","scaleup","pileup","piledown"]
	#~ suffixes = ["smeared"]
	lambdas = [10,16,22,28,34,40]
	interferences = ["Con","Des"]
	hels = ["LL","LR","RR"]

	#~ massBins = [1200,1400,1600,1800,2000,2200,2400,2600,2800,3000,3200,3400]
	massBins = [1200,1400,1600,1800,2000,2200,2400]
	signalYields = {}
	
	#~ names = ["default","resolution","scale","ID"]
	
	for label in labels:
		if "dimuon" in label:
			suffixes = suffixesMu
		else:
			suffixes = suffixesEle
		for suffix in suffixes:
			for histo in histos:
				for hel in hels:
					for interference in interferences:			
						model = interference+hel		
						if "dimuon" in label:
							name = "cito2mu"
						else:
							name = "cito2e"	
						fitFile = TFile("%s_%s_%s_inc_parametrization_fixinf.root"%(name,suffix,histo.lower()),"READ")
						for l in lambdas:
							if "dimuon" in label:
								name = "CITo2Mu_Lam%dTeV%s"%(l,model)
							else:	
								name = "CITo2E_Lam%dTeV%s"%(l,model)
							signalYields["%s_%s_%s"%(name,label,histo)] = {}
							for index, massBin in enumerate(massBins):
								function = fitFile.Get("fn_m%d_%s"%(massBin,model))
								fitR = fitFile.Get("fitR_m%d_%s"%(massBin,model))
								pars = fitR.GetParams()
								errs = fitR.Errors()
								function.SetParameter(0,pars[0])
								function.SetParameter(1,pars[1])
								function.SetParameter(2,pars[2])
								function.SetParError(0,errs[0])
								function.SetParError(1,errs[1])
								function.SetParError(2,errs[2])
								functionUnc = fitFile.Get("fn_unc_m%d_%s"%(massBin,model))
								print (function.Eval(l), function.Eval(100000), l, model, massBin)
								uncert = ((functionUnc.Eval(l)/function.Eval(l))**2 + (functionUnc.Eval(100000)/function.Eval(100000)))**0.5
								signalYields["%s_%s_%s"%(name,label,histo)][str(massBin)] = [(function.Eval(l)-function.Eval(100000)),uncert]
					if "dimuon" in label:
						name = "CITo2Mu_Lam%dTeV%s"%(l,hel)
						nameCon = "CITo2Mu_Lam%dTeV%s"%(l,"Con"+hel)
						nameDes = "CITo2Mu_Lam%dTeV%s"%(l,"Des"+hel)
					else:	
						name = "CITo2E_Lam%dTeV%s"%(l,hel)
						nameCon = "CITo2E_Lam%dTeV%s"%(l,"Con"+hel)
						nameDes = "CITo2E_Lam%dTeV%s"%(l,"Des"+hel)
					signalYields["%s_%s_%s"%(name,label,histo)] = {}	
					for index, massBin in enumerate(massBins):
						signalYields["%s_%s_%s"%(name,label,histo)][str(massBin)] = [(signalYields["%s_%s_%s"%(nameCon,label,histo)][str(massBin)][0] + signalYields["%s_%s_%s"%(nameDes,label,histo)][str(massBin)][0])/2,((signalYields["%s_%s_%s"%(nameCon,label,histo)][str(massBin)][1]**2 + signalYields["%s_%s_%s"%(nameDes,label,histo)][str(massBin)][1]**2))**0.5]
					
			



			print (signalYields)
			if "dimuon" in label:
				fileName = "signalYieldsSingleBin"
			else:
				fileName = "signalYieldsSingleBinEle"
			
			if suffix == "nominal":
				otherSuffix = "default"
			elif suffix == "scaledown":
				otherSuffix = "scaleDown"
			elif suffix == "scaleup":
				otherSuffix = "scaleUp"
			elif suffix == "smeared":
				otherSuffix = "resolution"
			elif suffix == "muonid":
				otherSuffix = "ID"
			elif suffix == "pileup":
				otherSuffix = "pileup"
			elif suffix == "piledown":
				otherSuffix = "piledown"
			else:
				print (suffix)
			outFilePkl = open("%s_%s.pkl"%(fileName,otherSuffix),"wb")
			pickle.dump(signalYields, outFilePkl)
			outFilePkl.close()		
	
	
			
							
main()
