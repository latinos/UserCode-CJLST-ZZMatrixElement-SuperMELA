#! /usr/bin/env python
import sys
import os
import re
import math
from scipy.special import erf
from ROOT import *
import ROOT
from array import array
from systematicsClassSMD import *
from inputReader import *

## ------------------------------------
##  card and workspace class
## ------------------------------------
##
##
## THIS VERSION IS FOR CARDS FOR HYPOTHESIS TEST WITH SUPERMELA COMBINE !!!
##
##
class datacardClass:

    def __init__(self):
    
        self.ID_4mu = 1
        self.ID_4e  = 2
        self.ID_2e2mu = 3    
        self.isFSR = True

    def loadIncludes(self):
        
        ROOT.gSystem.AddIncludePath("-I$ROOFITSYS/include/")
        ROOT.gSystem.AddIncludePath("-Iinclude/")
        ROOT.gROOT.ProcessLine(".L include/tdrstyle.cc")
        ROOT.gSystem.Load("libRooFit")
        ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
        ROOT.gSystem.Load("include/HiggsCSandWidth_cc.so")
        ROOT.gSystem.Load("include/HiggsCSandWidthSM4_cc.so")

    def reflectSystematics(self,nomShape,altShape):

        if(nomShape.GetNbinsX()!=altShape.GetNbinsX() or nomShape.GetNbinsY()!=altShape.GetNbinsY()):
            print "AHHHHHHHHHHH, templates don't have the same binning!!!!"
            return 0

        newAltShape = ROOT.TH2F(altShape)

        for x in range(1,nomShape.GetNbinsX()):
            for y in range(1,nomShape.GetNbinsY()):
                delta=altShape.GetBinContent(x,y)-nomShape.GetBinContent(x,y)
                newAltShape.SetBinContent(x,y,nomShape.GetBinContent(x,y)-delta)
                if(newAltShape.GetBinContent(x,y)<0.0):
                    newAltShape.SetBinContent(x,y,0.0)
            # done with loop over y bins
        #done with loop over x bins

        newAltShape.Scale(1.0/newAltShape.Integral())

        #check that no bins are zero

        for x in range(1,newAltShape.GetNbinsX()):
            for y in range(1,newAltShape.GetNbinsY()):
                if(newAltShape.GetBinContent(x,y)<0.000001):
                    newAltShape.SetBinContent(x,y,0.000001)
            # done with loop over y bins
        #done with loop over x bins

        newAltShape.Scale(1.0/newAltShape.Integral())

        return newAltShape

    # cross section filter for 7 TeV efficiency
    def csFilter(self,hmass):

        a = 80.85
        b = 50.42
        
        f = 0.5 + 0.5*erf( (hmass - a)/b )
        
        return f

    # cs x br function 
    def makeXsBrFunction(self,signalProc,rrvMH):
            
        procName = "ggH"
        if(signalProc == 0): procName = "ggH" #dummy, when you sum up all the 5 chans
        if(signalProc == 1): procName = "ggH"
        if(signalProc == 2): procName = "qqH"
        if(signalProc == 3): procName = "WH"
        if(signalProc == 4): procName = "ZH"
        if(signalProc == 5): procName = "ttH"

        
        
        
        channelName = ""
        if (self.channel == self.ID_4mu): channelName = "4mu"
        elif (self.channel == self.ID_4e): channelName = "4e"
        elif (self.channel == self.ID_2e2mu): channelName = "2e2mu"
        else: print "Input Error: Unknown channel! (4mu = 1, 4e = 2, 2e2mu = 3)" 

     
        
        myCSWrhf = HiggsCSandWidth()
        
        histXsBr = ROOT.TH1F("hsmxsbr_{0}_{1}".format(procName,channelName),"", 1781, 109.75, 1000)
        
        for i in range(1,1782):
            
            mHVal = histXsBr.GetBinCenter(i)
            BR = 0.0 
            if (self.channel == self.ID_2e2mu):
                BR = myCSWrhf.HiggsBR(13,mHVal)
            else:
                BR = myCSWrhf.HiggsBR(12,mHVal)

            if (signalProc==0):
                totXs=0
                for ch in range(1,6):
                    totXs+=myCSWrhf.HiggsCS(ch, mHVal, self.sqrts)
                histXsBr.SetBinContent(i, totXs * BR)
            else:
                histXsBr.SetBinContent(i, myCSWrhf.HiggsCS(signalProc, mHVal, self.sqrts) * BR)

            #print '\nmakeXsBrFunction : procName=',procName,'   signalProc=',signalProc,'  mH (input)=',rrvMH.getVal(),
            #print '   CS=',myCSWrhf.HiggsCS(signalProc, mHVal, self.sqrts),'   BR=',BR
            
        rdhname = "rdhXsBr_{0}_{1}_{2}".format(procName,self.channel,self.sqrts)
        rdhXsBr = RooDataHist(rdhname,rdhname, ROOT.RooArgList(rrvMH), histXsBr)  
        
        return rdhXsBr
    
    # return trueVar if testStatement else return falseVar
    def getVariable(self,trueVar,falseVar,testStatement):

        if (testStatement): 
            return trueVar
        else:
            return falseVar
    
    # main datacard and workspace function
    def makeCardsWorkspaces(self, theMH, theis2D, theOutputDir, theInputs,theTemplateDir="templates2D", theSpinHyp=-1, theIncludingError=False):

        ## --------------- SETTINGS AND DECLARATIONS --------------- ##
        DEBUG = False
        self.mH = theMH
        self.lumi = theInputs['lumi']
        self.sqrts = theInputs['sqrts']
        self.channel = theInputs['decayChannel']
        self.is2D = theis2D
        self.outputDir = theOutputDir
        self.sigMorph = theInputs['useCMS_zz4l_sigMELA']
        self.bkgMorph = theInputs['useCMS_zz4l_bkgMELA']
        self.templateDir = theTemplateDir
        self.spinHypID= theSpinHyp
	self.bIncludingError=theIncludingError
	self.spinHyp = theSpinHyp
        
        FactorizedShapes = False

        self.all_chan = theInputs['all']
        self.ggH_chan = theInputs['ggH']
        self.qqH_chan = theInputs['qqH']
        self.WH_chan = theInputs['WH']
        self.ZH_chan = theInputs['ZH']
        self.ttH_chan = theInputs['ttH']
        self.qqZZ_chan = theInputs['qqZZ']
        self.ggZZ_chan = theInputs['ggZZ']
        self.zjets_chan = theInputs['zjets']
        self.ttbar_chan = theInputs['ttbar']
        self.zbb_chan = theInputs['zbb']
        
        ## ---------------- SET PLOTTING STYLE ---------------- ## 
        ROOT.setTDRStyle(True)
        ROOT.gStyle.SetPalette(1)
        ROOT.gStyle.SetPadLeftMargin(0.16)        

        ## ---------------- VARIABLES FOR LATER --------------- ##
        self.bUseCBnoConvolution = False
        ForXSxBR = False

        myCSW = HiggsCSandWidth()
                
        ## ----------------- WIDTH AND RANGES ----------------- ##
        self.widthHVal =  myCSW.HiggsWidth(0,self.mH)
        if(self.widthHVal < 0.12):
            self.bUseCBnoConvolution = True
        self.isHighMass = False
        if self.mH >= 400:
            if theInputs['useHighMassReweightedShapes']:
                self.isHighMass = True
            else: print "useHighMassReweightedShapes set to FALSE, using non-reweighted shapes!"

            
        if(DEBUG): print "width: ",self.widthHVal
        
        self.windowVal = max( self.widthHVal, 1.0)
        lowside = 100.0
        if (self.mH >= 275):
            lowside = 180.0
        else:
            lowside = 100.0
        
        self.low_M = max( (self.mH - 20.*self.windowVal), lowside)
        self.high_M = min( (self.mH + 15.*self.windowVal), 1000)

        #self.low_M = 100.0
        #self.high_M = 800.0
       
        if (self.channel == self.ID_4mu): self.appendName = '4mu'
        elif (self.channel == self.ID_4e): self.appendName = '4e'
        elif (self.channel == self.ID_2e2mu): self.appendName = '2e2mu'
        else: print "Input Error: Unknown channel! (4mu = 1, 4e = 2, 2e2mu = 3)"

        self.isAltSig = False
        if (theInputs['doHypTest']):
            self.isAltSig = True
            
        if self.isAltSig and not self.all_chan :
            raise RuntimeError, "You asked to prepare DC and WS for Hyp Test but you did not want to sum over all signal channels. This is forbidden. Check inputs ! (it should have already send you this error message, strange that  you are here...)"

        if (self.isAltSig and not (self.is2D==1)):
            raise RunTimeError, "Cannot perform hypothesis testing without a 2D analysis, feature not supported yet. Exiting."
        

        self.appendHypType = theInputs['altHypLabel']
        if self.isAltSig and self.appendHypType=="" :
            self.appendHypType = "_ALT"
            
        
        ## ------------------------- SYSTEMATICS CLASSES ----------------------------- ##
    
        systematics = systematicsClass( self.mH, False, self.isFSR, theInputs)
        systematics_forXSxBR = systematicsClass( self.mH, True, self.isFSR,theInputs)

        ## -------------------------- SIGNAL SHAPE ----------------------------------- ##
    
        ##bins = 1000
        bins = 200
        

        CMS_zz4l_mass = ROOT.RooRealVar("CMS_zz4l_mass","CMS_zz4l_mass",self.low_M,self.high_M)
        CMS_zz4l_mass.setBins(bins,"fft") 

        superDiscVarName = "CMS_zz4l_smd"
        SD = ROOT.RooRealVar(superDiscVarName,superDiscVarName,0.0,1.0)

        self.LUMI = ROOT.RooRealVar("LUMI_{0:.0f}".format(self.sqrts),"LUMI_{0:.0f}".format(self.sqrts),self.lumi)
        self.LUMI.setConstant(True)
    
        self.MH = ROOT.RooRealVar("MH","MH",self.mH)
        self.MH.setConstant(True)
    
	# bIncludingError 
	CMS_zz4l_massErr = ROOT.RooRealVar("CMS_zz4l_massErr", "CMS_zz4l_massErr", 0.01*self.low_M/3., 0.01*self.high_M*5 );
	CMS_zz4l_massErr.setBins(100);

	# n2, alpha2 are right side parameters of DoubleCB
	# n, alpha are left side parameters of DoubleCB

        n_CB_d = 0.0
        alpha_CB_d = 0.0
        n2_CB_d = 0.0
        alpha2_CB_d = 0.0
        mean_CB_d = 0.0
        sigma_CB_d = 0.0
        mean_BW_d = self.mH
        gamma_BW_d = 0.0
        
        if(self.all_chan):
            rdhXsBrFuncV_1 = self.makeXsBrFunction(0,self.MH)
        else:
            rdhXsBrFuncV_1 = self.makeXsBrFunction(1,self.MH)
        rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}".format("ggH",self.channel,self.sqrts)
        rhfXsBrFuncV_1 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), rdhXsBrFuncV_1, 1)
        
        rdhXsBrFuncV_2 = self.makeXsBrFunction(2,self.MH)
        rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}".format("VBF",self.channel,self.sqrts)
        rhfXsBrFuncV_2 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), rdhXsBrFuncV_2, 1)
        
        rdhXsBrFuncV_3 = self.makeXsBrFunction(3,self.MH)
        rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}".format("WH",self.channel,self.sqrts)
        rhfXsBrFuncV_3 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), rdhXsBrFuncV_3, 1)
        
        rdhXsBrFuncV_4 = self.makeXsBrFunction(4,self.MH)
        rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}".format("ZH",self.channel,self.sqrts)
        rhfXsBrFuncV_4 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), rdhXsBrFuncV_4, 1)
        
        rdhXsBrFuncV_5 = self.makeXsBrFunction(5,self.MH)
        rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}".format("ttH",self.channel,self.sqrts)
        rhfXsBrFuncV_5 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), rdhXsBrFuncV_5, 1)

    
        ## -------- Variable Definitions -------- ##
    
        name = "CMS_zz4l_mean_e_sig"
        CMS_zz4l_mean_e_sig = ROOT.RooRealVar(name,"CMS_zz4l_mean_e_sig",0.0,-10.0,10.0)
        name = "CMS_zz4l_sigma_e_sig"
        CMS_zz4l_sigma_e_sig = ROOT.RooRealVar(name,"CMS_zz4l_sigma_sig",3.0,0.0,30.0)
        
        name = "CMS_zz4l_mean_m_sig"
        CMS_zz4l_mean_m_sig = ROOT.RooRealVar(name,"CMS_zz4l_mean_sig",0.0,-10.0,10.0)
        name = "CMS_zz4l_sigma_m_sig"
        CMS_zz4l_sigma_m_sig = ROOT.RooRealVar(name,"CMS_zz4l_sigma_sig",3.0,0.0,30.0)
        
        name = "CMS_zz4l_alpha2_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_alpha2 = ROOT.RooRealVar(name,"CMS_zz4l_alpha2",1.,-10.,10.)
        name = "CMS_zz4l_n2_sig_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_n2 = ROOT.RooRealVar(name,"CMS_zz4l_n2",2.,-10.,10.)
        name = "CMS_zz4l_alpha_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_alpha = ROOT.RooRealVar(name,"CMS_zz4l_alpha",1.,-10.,10.)
        name = "CMS_zz4l_n_sig_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_n = ROOT.RooRealVar(name,"CMS_zz4l_n",2.,-10.,10.)
        name = "CMS_zz4l_mean_BW_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_mean_BW = ROOT.RooRealVar(name,"CMS_zz4l_mean_BW",self.mH,self.low_M,self.high_M)
        name = "CMS_zz4l_gamma_sig_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_gamma = ROOT.RooRealVar(name,"CMS_zz4l_gamma",10.,0.001,1000.)
        name = "CMS_zz4l_widthScale_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_widthScale = ROOT.RooRealVar(name,"CMS_zz4l_widthScale",1.0)
        one = ROOT.RooRealVar("one","one",1.0)
        one.setConstant(True)
    
        CMS_zz4l_mean_BW.setVal( mean_BW_d )
        CMS_zz4l_gamma.setVal(0)
        CMS_zz4l_mean_e_sig.setVal(0)
        CMS_zz4l_sigma_e_sig.setVal(0)
        CMS_zz4l_mean_m_sig.setVal(0)
        CMS_zz4l_sigma_m_sig.setVal(0)
        CMS_zz4l_alpha.setVal(0)
        CMS_zz4l_n.setVal(0)
        CMS_zz4l_alpha2.setVal(0)
        CMS_zz4l_n2.setVal(0)
    
        CMS_zz4l_widthScale.setConstant(True)
        #CMS_zz4l_alpha.setConstant(True)  # also read from input file
        CMS_zz4l_mean_BW.setConstant(True)
        #CMS_zz4l_gamma_BW.setConstant(True)

        print "HEEERRRRRRRRRRRRRRRRREEEEEEE"

        print "mean_BW ", CMS_zz4l_mean_BW.getVal()
        print "gamma_BW ", CMS_zz4l_gamma.getVal()
        print "mean_e_sig ", CMS_zz4l_mean_e_sig.getVal()
        print "sigma_e ", CMS_zz4l_sigma_e_sig.getVal()
        print "mean_m_sig ",CMS_zz4l_mean_m_sig.getVal()
        print "sigma_m ", CMS_zz4l_sigma_m_sig.getVal()
        print "alpha ", CMS_zz4l_alpha.getVal()
        print "n ", CMS_zz4l_n.getVal()
        print "alpha2 ", CMS_zz4l_alpha2.getVal()
        print "n2 ", CMS_zz4l_n2.getVal()

                                                                


        ## -------------------- RooFormulaVar's -------------------- ##
        rfv_n_CB = ROOT.RooFormulaVar()
        rfv_alpha_CB = ROOT.RooFormulaVar()
        rfv_n2_CB = ROOT.RooFormulaVar()
        rfv_alpha2_CB = ROOT.RooFormulaVar()
        rfv_mean_CB = ROOT.RooFormulaVar()
        rfv_sigma_CB = ROOT.RooFormulaVar()
        
        name = "CMS_zz4l_n_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
        rfv_n_CB = ROOT.RooFormulaVar(name,"("+theInputs['n_CB_shape']+")"+"*(1+@1)",ROOT.RooArgList(self.MH,CMS_zz4l_n))

        name = "CMS_zz4l_alpha_{0:.0f}_centralValue".format(self.channel)
        rfv_alpha_CB = ROOT.RooFormulaVar(name,theInputs['alpha_CB_shape'], ROOT.RooArgList(self.MH))

        name = "CMS_zz4l_n2_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
        rfv_n2_CB = ROOT.RooFormulaVar(name,"("+theInputs['n2_CB_shape']+")",ROOT.RooArgList(self.MH))

        name = "CMS_zz4l_alpha2_{0:.0f}_centralValue".format(self.channel)
        rfv_alpha2_CB = ROOT.RooFormulaVar(name,theInputs['alpha2_CB_shape'], ROOT.RooArgList(self.MH))

        if (DEBUG): print " DEBUG *********  ", theInputs['mean_CB_shape']
        name = "CMS_zz4l_mean_sig_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
        if (DEBUG): print " DEBUG *********  Name of rfv_mean_CB: ",name
        if (self.channel == self.ID_4mu) :
            rfv_mean_CB = ROOT.RooFormulaVar(name,"("+theInputs['mean_CB_shape']+")"+"+@0*@1", ROOT.RooArgList(self.MH, CMS_zz4l_mean_m_sig))
        elif (self.channel == self.ID_4e) :
            rfv_mean_CB = ROOT.RooFormulaVar(name,"("+theInputs['mean_CB_shape']+")"+"+@0*@1", ROOT.RooArgList(self.MH, CMS_zz4l_mean_e_sig))
        elif (self.channel == self.ID_2e2mu) :
            rfv_mean_CB = ROOT.RooFormulaVar(name,"("+theInputs['mean_CB_shape']+")"+"+ @0*@1 + @0*@2", ROOT.RooArgList(self.MH, CMS_zz4l_mean_m_sig,CMS_zz4l_mean_e_sig))
        if (DEBUG): print " DEBUG *********  Valuee of rfv_mean_CB: ",rfv_mean_CB.getVal()

        name = "CMS_zz4l_sigma_sig_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
        if (self.channel == self.ID_4mu) :
            rfv_sigma_CB = ROOT.RooFormulaVar(name,"("+theInputs['sigma_CB_shape']+")"+"*(1+@1)", ROOT.RooArgList(self.MH, CMS_zz4l_sigma_m_sig))
        elif (self.channel == self.ID_4e) :
            rfv_sigma_CB = ROOT.RooFormulaVar(name,"("+theInputs['sigma_CB_shape']+")"+"*(1+@1)", ROOT.RooArgList(self.MH, CMS_zz4l_sigma_e_sig))
        elif (self.channel == self.ID_2e2mu) :
            rfv_sigma_CB = ROOT.RooFormulaVar(name,"("+theInputs['sigma_CB_shape']+")"+"*(1+@1)*(1+@2)", ROOT.RooArgList(self.MH, CMS_zz4l_sigma_m_sig,CMS_zz4l_sigma_e_sig))


        name = "CMS_zz4l_gamma_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
        rfv_gamma_BW = ROOT.RooFormulaVar(name,"("+theInputs['gamma_BW_shape_HM']+")"+"*(1+@1)",ROOT.RooArgList(self.MH,CMS_zz4l_gamma))

        if (DEBUG): print " DEBUG *********  ", theInputs['sigma_CB_shape'] 

        print "n_CB ", rfv_n_CB.getVal()
        print "alpha_CB ", rfv_alpha_CB.getVal()
        print "n2_CB ", rfv_n2_CB.getVal()
        print "alpha2_CB ", rfv_alpha2_CB.getVal()
        print "mean_CB ", rfv_mean_CB.getVal()
        print "sigma_CB ", rfv_sigma_CB.getVal()
        print "gamma_BW ", rfv_gamma_BW.getVal()    

        
        CMS_zz4l_mean_sig_NoConv = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_NoConv_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts),"@0+@1", ROOT.RooArgList(rfv_mean_CB, self.MH))

        print "mean_sig_NoConv ", CMS_zz4l_mean_sig_NoConv.getVal()
        
        ## --------------------- SHAPE FUNCTIONS ---------------------- ##
        print 'Value of CMS_zz4l_mass ',CMS_zz4l_mass.getVal()
        signalCB_ggH = ROOT.RooDoubleCB("signalCB_ggH","signalCB_ggH",CMS_zz4l_mass, self.getVariable(CMS_zz4l_mean_sig_NoConv,rfv_mean_CB, self.bUseCBnoConvolution) , self.getVariable(CMS_zz4l_massErr,rfv_sigma_CB, self.bIncludingError),rfv_alpha_CB,rfv_n_CB, rfv_alpha2_CB, rfv_n2_CB)
        print 'Value of signalCB_ggH', signalCB_ggH.getVal()
        #Low mass pdf
        signalBW_ggH = ROOT.RooRelBWUFParam("signalBW_ggH", "signalBW_ggH",CMS_zz4l_mass,CMS_zz4l_mean_BW,CMS_zz4l_widthScale)
        print 'Value of signalBW_ggH', signalBW_ggH.getVal()

        sig_ggH =  ROOT.RooFFTConvPdf("sig_ggH","BW (X) CB",CMS_zz4l_mass,signalBW_ggH,signalCB_ggH, 2)
  
        
        ## Buffer fraction for cyclical behavior
        sig_ggH.setBufferFraction(0.2)
        print 'Value of sig_ggH', sig_ggH.getVal()
        ## --------------------------- SuperMELA 2D PDFS ------------------------- ##


        discVarName = "CMS_zz4l_pseudoKD"
        if(self.spinHyp == 0):
            discVarName = "CMS_zz4l_pseudoKD"
        elif (self.spinHyp == 1):
            discVarName = "CMS_zz4l_graviKD"
        elif (self.spinHyp == 2):
            discVarName = "CMS_zz4l_p0hplusKD"
        elif (self.spinHyp == 3):
            discVarName = "CMS_zz4l_p1plusKD"
        elif (self.spinHyp == 4):
            discVarName = "CMS_zz4l_p1minusKD"
        elif (self.spinHyp == 5):
            discVarName = "CMS_zz4l_qqgraviKD"
        else :
            discVarName = "CMS_zz4l_pseudoKD"

        print '++++ SuperMELA 2D PDFS using discriminat named :',discVarName,'  +++++'
        D = ROOT.RooRealVar(discVarName,discVarName,0,1)
    
        templateSigName = "{0}/Dsignal_{1}.root".format(self.templateDir ,self.appendName)
        
        sigTempFile = ROOT.TFile(templateSigName)
        sigTemplate = sigTempFile.Get("h_superDpsD")
        sigTemplate_syst1Up = sigTempFile.Get("h_superDpsD_LeptScaleUp")
        sigTemplate_syst1Down = sigTempFile.Get("h_superDpsD_LeptScaleDown")
        sigTemplate_syst2Up = sigTempFile.Get("h_superDpsD_LeptSmearUp")
        sigTemplate_syst2Down = sigTempFile.Get("h_superDpsD_LeptSmearDown")
            
        TemplateName = "sigTempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(SD,D),sigTemplate)
        TemplateName = "sigTempDataHist_syst1Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTempDataHist_syst1Up = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(SD,D),sigTemplate_syst1Up)
        TemplateName = "sigTempDataHist_syst1Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTempDataHist_syst1Down =ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(SD,D),sigTemplate_syst1Down)
        TemplateName = "sigTempDataHist_syst2Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTempDataHist_syst2Up = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(SD,D),sigTemplate_syst1Up)
        TemplateName = "sigTempDataHist_syst2Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTempDataHist_syst2Down =ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(SD,D),sigTemplate_syst1Down)

        TemplateName = "sigTemplatePdf_ggH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplatePdf_ggH = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(SD,D),sigTempDataHist)
        TemplateName = "sigTemplatePdf_ggH_syst1Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplatePdf_ggH_syst1Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(SD,D),sigTempDataHist_syst1Up)
        TemplateName = "sigTemplatePdf_ggH_syst1Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplatePdf_ggH_syst1Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(SD,D),sigTempDataHist_syst1Down)
        TemplateName = "sigTemplatePdf_ggH_syst2Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplatePdf_ggH_syst2Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(SD,D),sigTempDataHist_syst2Up)
        TemplateName = "sigTemplatePdf_ggH_syst2Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplatePdf_ggH_syst2Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(SD,D),sigTempDataHist_syst2Down)
           

        if(self.isAltSig):
            #only ggH because if we do hypothesis testing we sum up over the channels in any case                                           
            templateSigName = "{0}/Dsignal{2}_{1}.root".format(self.templateDir,self.appendName, self.appendHypType)
            print 'Taking 2D template for ALT signal from ',templateSigName
            sigTempFile = ROOT.TFile(templateSigName)
            sigTemplate = sigTempFile.Get("h_superDpsD")
            #systematics                                                                                                                  
            sigTemplate_syst1Up = sigTempFile.Get("h_superDpsD_LeptScaleUp")
            sigTemplate_syst1Down = sigTempFile.Get("h_superDpsD_LeptScaleDown")
            sigTemplate_syst2Up = sigTempFile.Get("h_superDpsD_LeptSmearUp")
            sigTemplate_syst2Down = sigTempFile.Get("h_superDpsD_LeptSmearDown")
            
            TemplateName = "sigTempDataHist_{0:.0f}_{1:.0f}{2}".format(self.channel,self.sqrts, self.appendHypType)
            sigTempDataHist_ALT = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(SD,D),sigTemplate)
            TemplateName = "sigTempDataHist_syst1Up_{0:.0f}_{1:.0f}{2}".format(self.channel,self.sqrts, self.appendHypType)
            sigTempDataHist_ALT_syst1Up = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(SD,D),sigTemplate_syst1Up)
            TemplateName = "sigTempDataHist_syst1Down_{0:.0f}_{1:.0f}{2}".format(self.channel,self.sqrts, self.appendHypType)
            sigTempDataHist_ALT_syst1Down = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(SD,D),sigTemplate_syst1Down)
            TemplateName = "sigTempDataHist_syst2Up_{0:.0f}_{1:.0f}{2}".format(self.channel,self.sqrts, self.appendHypType)
            sigTempDataHist_ALT_syst2Up = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(SD,D),sigTemplate_syst2Up)
            TemplateName = "sigTempDataHist_syst2Down_{0:.0f}_{1:.0f}{2}".format(self.channel,self.sqrts, self.appendHypType)
            sigTempDataHist_ALT_syst2Down = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(SD,D),sigTemplate_syst2Down)
            
            
            TemplateName = "sigTemplatePdf_ggH{2}_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts, self.appendHypType)
            sigTemplatePdf_ggH_ALT = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(SD,D),sigTempDataHist_ALT)
            TemplateName = "sigTemplatePdf_ggH{2}_syst1Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts, self.appendHypType)
            sigTemplatePdf_ggH_ALT_syst1Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(SD,D),sigTempDataHist_ALT_syst1Up)
            TemplateName = "sigTemplatePdf_ggH{2}_syst1Down_{0:.0f}_{1:.0f}{2}".format(self.channel,self.sqrts, self.appendHypType)
            sigTemplatePdf_ggH_ALT_syst1Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(SD,D),sigTempDataHist_ALT_syst1Down)
            TemplateName = "sigTemplatePdf_ggH{2}_syst2Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts, self.appendHypType)
            sigTemplatePdf_ggH_ALT_syst2Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(SD,D),sigTempDataHist_ALT_syst2Up)
            TemplateName = "sigTemplatePdf_ggH{2}_syst2Down_{0:.0f}_{1:.0f}{2}".format(self.channel,self.sqrts, self.appendHypType)
            sigTemplatePdf_ggH_ALT_syst2Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(SD,D),sigTempDataHist_ALT_syst2Down)

             
            ###Shape systematics for signal
        funcList1_ggH = ROOT.RooArgList()  
        funcList1_ggH_ALT = ROOT.RooArgList()
       
        if(self.sigMorph): #this switches on the inclusion of shape systematics
            
            funcList1_ggH.add(sigTemplatePdf_ggH)
            funcList1_ggH.add(sigTemplatePdf_ggH_syst1Up)
            funcList1_ggH.add(sigTemplatePdf_ggH_syst1Down)  

            funcList1_ggH.add(sigTemplatePdf_ggH_syst2Up)
            funcList1_ggH.add(sigTemplatePdf_ggH_syst2Down)  

            if(self.isAltSig):
                funcList1_ggH_ALT.add(sigTemplatePdf_ggH_ALT)
                funcList1_ggH_ALT.add(sigTemplatePdf_ggH_ALT_syst1Up)
                funcList1_ggH_ALT.add(sigTemplatePdf_ggH_ALT_syst1Down)
                funcList1_ggH_ALT.add(sigTemplatePdf_ggH_ALT_syst2Up)
                funcList1_ggH_ALT.add(sigTemplatePdf_ggH_ALT_syst2Down)
        else:
            
            funcList1_ggH.add(sigTemplatePdf_ggH)
            if(self.isAltSig):
                funcList1_ggH_ALT.add(sigTemplatePdf_ggH_ALT)
    
        morphSigVarName = "CMS_zz4l_smd_leptScale_sig_{0:.0f}".format(self.channel)
        syst1MorphSig = ROOT.RooRealVar(morphSigVarName,morphSigVarName,0,-20,20)
        morphSigVarName = "CMS_zz4l_smd_leptResol_sig_{0:.0f}".format(self.channel)
        syst2MorphSig = ROOT.RooRealVar(morphSigVarName,morphSigVarName,0,-20,20)
        if(self.sigMorph):
            syst1MorphSig.setConstant(False)
            syst2MorphSig.setConstant(False)
        else:
            syst1MorphSig.setConstant(True)
            syst2MorphSig.setConstant(True)
        
        morphVarListSig1 = ROOT.RooArgList()
###        morphVarListSig2 = ROOT.RooArgList() ## just one morphing for all signal processes (fully correlated systs) 
    
        if(self.sigMorph):
            morphVarListSig1.add(syst1MorphSig)  ## just one morphing for all signal processes (fully correlated systs)
###            morphVarListSig2.add(syst2MorphSig)
            morphVarListSig1.add(syst2MorphSig)
            
        TemplateName = "sigTemplateMorphPdf_ggH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplateMorphPdf_ggH = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,SD,D,false,funcList1_ggH,morphVarListSig1,1.0,1)
        
        if(self.isAltSig):
            TemplateName = "sigTemplateMorphPdf_ggH{2}_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts,self.appendHypType)
            sigTemplateMorphPdf_ggH_ALT = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,SD,D,false,funcList1_ggH_ALT,morphVarListSig1,1.0,1)


        sigCB2d_ggH =signalCB_ggH.Clone("sigCB2d_ggH")

        if(self.isAltSig):
            sigCB2d_ggH_ALT =signalCB_ggH.Clone("sigCB2d_ggH{0}".format(self.appendHypType))



        ## --------------------------- superMELA 1D PDFS ------------------------- ##     
    
        templateSDSigName = "{0}/Dsignal_{1}.root".format(self.templateDir ,self.appendName)
        sigTempSDFile = ROOT.TFile(templateSDSigName)
        sigTemplateSD = sigTempSDFile.Get("h_superD")
###        sigTemplateSD = sigTempSDFile.Get("h_superD_mod_rndm")
###        sigTemplateSD = sigTempSDFile.Get("h_superDfromProjX") 
        print 'SuperMELA 1D template has mean = ',sigTemplateSD.GetMean()
        
        TemplateSDName = "sigTempSDDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTempSDDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(SD),sigTemplateSD)
        
        TemplateSDName = "sigTemplateSDPdf_ggH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplateSDPdf_ggH = ROOT.RooHistPdf(TemplateSDName,TemplateSDName,ROOT.RooArgSet(SD),sigTempSDDataHist)
        
        TemplateSDName = "sigTemplateSDPdf_VBF_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplateSDPdf_VBF = ROOT.RooHistPdf(TemplateSDName,TemplateSDName,RooArgSet(SD),sigTempSDDataHist)
        
        TemplateSDName = "sigTemplateSDPdf_WH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplateSDPdf_WH = ROOT.RooHistPdf(TemplateSDName,TemplateSDName,ROOT.RooArgSet(SD),sigTempSDDataHist)
        
        TemplateSDName = "sigTemplateSDPdf_ZH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplateSDPdf_ZH = ROOT.RooHistPdf(TemplateSDName,TemplateSDName,ROOT.RooArgSet(SD),sigTempSDDataHist)
        
        TemplateSDName = "sigTemplateSDPdf_ZH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplateSDPdf_ttH = ROOT.RooHistPdf(TemplateSDName,TemplateSDName,ROOT.RooArgSet(SD),sigTempSDDataHist)
        print sigTemplateSDPdf_ttH

        ##--------------##

        ## -------------------------- BACKGROUND SHAPES ---------------------------------- ##
        print 'BACKGROUND SHAPES'
        ## qqZZ contribution
        name = "CMS_qqzzbkg_a0_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a0 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a0",115.3,0.,200.)
        name = "CMS_qqzzbkg_a1_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a1 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a1",21.96,0.,200.)
        name = "CMS_qqzzbkg_a2_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a2 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a2",122.8,0.,200.)
        name = "CMS_qqzzbkg_a3_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a3 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a3",0.03479,0.,1.)
        name = "CMS_qqzzbkg_a4_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a4 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a4",185.5,0.,200.)
        name = "CMS_qqzzbkg_a5_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a5 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a5",12.67,0.,200.)
        name = "CMS_qqzzbkg_a6_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a6 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a6",34.81,0.,100.)
        name = "CMS_qqzzbkg_a7_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a7 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a7",0.1393,0.,1.)
        name = "CMS_qqzzbkg_a8_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a8 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a8",66.,0.,200.)
        name = "CMS_qqzzbkg_a9_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a9 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a9",0.07191,0.,1.)
        name = "CMS_qqzzbkg_a10_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts )
        CMS_qqzzbkg_a10 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a10",94.11,0.,200.)
        name = "CMS_qqzzbkg_a11_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts )
        CMS_qqzzbkg_a11 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a11",-5.111,-100.,100.)
        name = "CMS_qqzzbkg_a12_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts )
        CMS_qqzzbkg_a12 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a12",4834,0.,10000.)
        name = "CMS_qqzzbkg_a13_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts )
        CMS_qqzzbkg_a13 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a13",0.2543,0.,1.)
        

        if (DEBUG) :
            print "qqZZshape_a0 = ",theInputs['qqZZshape_a0']
            print "qqZZshape_a1 = ",theInputs['qqZZshape_a1']
            print "qqZZshape_a2 = ",theInputs['qqZZshape_a2']
            print "qqZZshape_a3 = ",theInputs['qqZZshape_a3']
            print "qqZZshape_a4 = ",theInputs['qqZZshape_a4']
            print "qqZZshape_a5 = ",theInputs['qqZZshape_a5']
            print "qqZZshape_a6 = ",theInputs['qqZZshape_a6']
            print "qqZZshape_a7 = ",theInputs['qqZZshape_a7']
            print "qqZZshape_a8 = ",theInputs['qqZZshape_a8']
            print "qqZZshape_a9 = ",theInputs['qqZZshape_a9']
            print "qqZZshape_a10 = ",theInputs['qqZZshape_a10']
            print "qqZZshape_a11 = ",theInputs['qqZZshape_a11']
            print "qqZZshape_a12 = ",theInputs['qqZZshape_a12']
            print "qqZZshape_a13 = ",theInputs['qqZZshape_a13']

        
        CMS_qqzzbkg_a0.setVal(theInputs['qqZZshape_a0'])
        CMS_qqzzbkg_a1.setVal(theInputs['qqZZshape_a1'])
        CMS_qqzzbkg_a2.setVal(theInputs['qqZZshape_a2'])
        CMS_qqzzbkg_a3.setVal(theInputs['qqZZshape_a3'])
        CMS_qqzzbkg_a4.setVal(theInputs['qqZZshape_a4'])
        CMS_qqzzbkg_a5.setVal(theInputs['qqZZshape_a5'])
        CMS_qqzzbkg_a6.setVal(theInputs['qqZZshape_a6'])
        CMS_qqzzbkg_a7.setVal(theInputs['qqZZshape_a7'])
        CMS_qqzzbkg_a8.setVal(theInputs['qqZZshape_a8'])
        CMS_qqzzbkg_a9.setVal(theInputs['qqZZshape_a9'])
        CMS_qqzzbkg_a10.setVal(theInputs['qqZZshape_a10'])
        CMS_qqzzbkg_a11.setVal(theInputs['qqZZshape_a11'])
        CMS_qqzzbkg_a12.setVal(theInputs['qqZZshape_a12'])
        CMS_qqzzbkg_a13.setVal(theInputs['qqZZshape_a13'])
        
        CMS_qqzzbkg_a0.setConstant(True)
        CMS_qqzzbkg_a1.setConstant(True)
        CMS_qqzzbkg_a2.setConstant(True)
        CMS_qqzzbkg_a3.setConstant(True)
        CMS_qqzzbkg_a4.setConstant(True)
        CMS_qqzzbkg_a5.setConstant(True)
        CMS_qqzzbkg_a6.setConstant(True)
        CMS_qqzzbkg_a7.setConstant(True)
        CMS_qqzzbkg_a8.setConstant(True)
        CMS_qqzzbkg_a9.setConstant(True)
        CMS_qqzzbkg_a10.setConstant(True)
        CMS_qqzzbkg_a11.setConstant(True)
        CMS_qqzzbkg_a12.setConstant(True)
        CMS_qqzzbkg_a13.setConstant(True)
        
        bkg_qqzz = ROOT.RooqqZZPdf_v2("bkg_qqzzTmp","bkg_qqzzTmp",CMS_zz4l_mass,CMS_qqzzbkg_a0,CMS_qqzzbkg_a1,CMS_qqzzbkg_a2,CMS_qqzzbkg_a3,CMS_qqzzbkg_a4,CMS_qqzzbkg_a5,CMS_qqzzbkg_a6,CMS_qqzzbkg_a7,CMS_qqzzbkg_a8,CMS_qqzzbkg_a9,CMS_qqzzbkg_a10,CMS_qqzzbkg_a11,CMS_qqzzbkg_a12,CMS_qqzzbkg_a13)
        
        ## ggZZ contribution
        name = "CMS_ggzzbkg_a0_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        CMS_ggzzbkg_a0 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a0",115.3,0.,200.)
        name = "CMS_ggzzbkg_a1_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        CMS_ggzzbkg_a1 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a1",21.96,0.,200.)
        name = "CMS_ggzzbkg_a2_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        CMS_ggzzbkg_a2 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a2",122.8,0.,200.)
        name = "CMS_ggzzbkg_a3_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        CMS_ggzzbkg_a3 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a3",0.03479,0.,1.)
        name = "CMS_ggzzbkg_a4_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts )
        CMS_ggzzbkg_a4 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a4",185.5,0.,200.)
        name = "CMS_ggzzbkg_a5_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        CMS_ggzzbkg_a5 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a5",12.67,0.,200.)
        name = "CMS_ggzzbkg_a6_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        CMS_ggzzbkg_a6 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a6",34.81,0.,100.)
        name = "CMS_ggzzbkg_a7_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        CMS_ggzzbkg_a7 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a7",0.1393,0.,1.)
        name = "CMS_ggzzbkg_a8_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        CMS_ggzzbkg_a8 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a8",66.,0.,200.)
        name = "CMS_ggzzbkg_a9_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts )
        CMS_ggzzbkg_a9 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a9",0.07191,0.,1.)
        
        CMS_ggzzbkg_a0.setVal(theInputs['ggZZshape_a0'])
        CMS_ggzzbkg_a1.setVal(theInputs['ggZZshape_a1'])
        CMS_ggzzbkg_a2.setVal(theInputs['ggZZshape_a2'])
        CMS_ggzzbkg_a3.setVal(theInputs['ggZZshape_a3'])
        CMS_ggzzbkg_a4.setVal(theInputs['ggZZshape_a4'])
        CMS_ggzzbkg_a5.setVal(theInputs['ggZZshape_a5'])
        CMS_ggzzbkg_a6.setVal(theInputs['ggZZshape_a6'])
        CMS_ggzzbkg_a7.setVal(theInputs['ggZZshape_a7'])
        CMS_ggzzbkg_a8.setVal(theInputs['ggZZshape_a8'])
        CMS_ggzzbkg_a9.setVal(theInputs['ggZZshape_a9'])
        
        CMS_ggzzbkg_a0.setConstant(True)
        CMS_ggzzbkg_a1.setConstant(True)
        CMS_ggzzbkg_a2.setConstant(True)
        CMS_ggzzbkg_a3.setConstant(True)
        CMS_ggzzbkg_a4.setConstant(True)
        CMS_ggzzbkg_a5.setConstant(True)
        CMS_ggzzbkg_a6.setConstant(True)
        CMS_ggzzbkg_a7.setConstant(True)
        CMS_ggzzbkg_a8.setConstant(True)
        CMS_ggzzbkg_a9.setConstant(True)

        if (DEBUG) :
            print "ggZZshape_a0 = ",theInputs['ggZZshape_a0']
            print "ggZZshape_a1 = ",theInputs['ggZZshape_a1']
            print "ggZZshape_a2 = ",theInputs['ggZZshape_a2']
            print "ggZZshape_a3 = ",theInputs['ggZZshape_a3']
            print "ggZZshape_a4 = ",theInputs['ggZZshape_a4']
            print "ggZZshape_a5 = ",theInputs['ggZZshape_a5']
            print "ggZZshape_a6 = ",theInputs['ggZZshape_a6']
            print "ggZZshape_a7 = ",theInputs['ggZZshape_a7']
            print "ggZZshape_a8 = ",theInputs['ggZZshape_a8']
            print "ggZZshape_a9 = ",theInputs['ggZZshape_a9']
                   
        
        bkg_ggzz = ROOT.RooggZZPdf_v2("bkg_ggzzTmp","bkg_ggzzTmp",CMS_zz4l_mass,CMS_ggzzbkg_a0,CMS_ggzzbkg_a1,CMS_ggzzbkg_a2,CMS_ggzzbkg_a3,CMS_ggzzbkg_a4,CMS_ggzzbkg_a5,CMS_ggzzbkg_a6,CMS_ggzzbkg_a7,CMS_ggzzbkg_a8,CMS_ggzzbkg_a9)
    
        ## Reducible backgrounds
        val_meanL = float(theInputs['zjetsShape_mean'])
        val_sigmaL = float(theInputs['zjetsShape_sigma'])
 
        name = "mlZjet_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        mlZjet = ROOT.RooRealVar(name,"mean landau Zjet",val_meanL)
        name = "slZjet_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        slZjet = ROOT.RooRealVar(name,"sigma landau Zjet",val_sigmaL)
        bkg_zjets = ROOT.RooLandau("bkg_zjetsTmp","bkg_zjetsTmp",CMS_zz4l_mass,mlZjet,slZjet) 



      ## ----------------- SuperMELA 2D BACKGROUND SHAPES --------------- ##
        print 'SUPERMELA 2D BACKGROUND SHAPES'
        templateBkgName = "{0}/Dbackground_qqZZ_{1}.root".format(self.templateDir ,self.appendName)
        bkgTempFile = ROOT.TFile(templateBkgName)
        bkgTemplate = bkgTempFile.Get("h_superDpsD")
        
        TemplateName = "bkgTempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        bkgTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(SD,D),bkgTemplate)
        TemplateName = "zjetsTempDataHist_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        bkgTempDataHist_zjetsUp = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(SD,D),bkgTemplate)

        templateggBkgName = "{0}/Dbackground_ggZZ_{1}.root".format(self.templateDir ,self.appendName)
        ggbkgTempFile = ROOT.TFile(templateggBkgName)
        ggbkgTemplate = ggbkgTempFile.Get("h_superDpsD")
        TemplateName = "ggbkgTempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        ggbkgTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(SD,D),ggbkgTemplate)

        TemplateName = "bkgTemplatePdf_qqzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        bkgTemplatePdf_qqzz = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(SD,D),bkgTempDataHist)
        TemplateName = "bkgTemplatePdf_ggzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        
        bkgTemplatePdf_ggzz = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(SD,D),ggbkgTempDataHist)

        templateBkgName = "{0}/Dbackground_ZJetsCR_AllChans.root".format(self.templateDir)
        zjetsTempFile = ROOT.TFile(templateBkgName)
        zjetsTemplate = zjetsTempFile.Get("h_superDpsD")
        TemplateName = "zjetsTempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        zjetsTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(SD,D),zjetsTemplate)



        TemplateName = "bkgTemplatePdf_zjets_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        bkgTemplatePdf_zjets = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(SD,D),bkgTempDataHist)
#        bkgTemplatePdf_zjets = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(SD,D),zjetsTempDataHist)
        TemplateName = "bkgTemplatePdf_zjets_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        bkgTemplatePdf_zjets_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(SD,D),bkgTempDataHist_zjetsUp)


        zjetsTemplateDown = self.reflectSystematics(bkgTemplate,zjetsTemplate)
        #zjetsTemplateDown.setName("zjetsTemplateDown_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts))
        #zjetsTemplateDown.setTitle("zjetsTemplateDown_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts))
        TemplateName = "zjetsTempDataHistDown_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        zjetsTempDataHistDown = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(SD,D),zjetsTemplateDown)
        TemplateName = "bkgTemplatePdf_zjets_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        bkgTemplatePdf_zjets_Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(SD,D),bkgTempDataHist)
                

        funcList_zjets = ROOT.RooArgList()  
        morphBkgVarName =  "CMS_zz4l_smd_zjets_bkg_{0:.0f}".format(self.channel)    
        alphaMorphBkg = ROOT.RooRealVar(morphBkgVarName,morphBkgVarName,0,-20,20)
        morphVarListBkg = ROOT.RooArgList()
        
        if(self.bkgMorph):
            funcList_zjets.add(bkgTemplatePdf_zjets)
            funcList_zjets.add(bkgTemplatePdf_zjets_Up)
            funcList_zjets.add(bkgTemplatePdf_zjets_Down)  
            alphaMorphBkg.setConstant(False)
            morphVarListBkg.add(alphaMorphBkg)  
        else:
            funcList_zjets.add(bkgTemplatePdf_zjets)
            alphaMorphBkg.setConstant(True)


        TemplateName = "bkgTemplateMorphPdf_qqzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        bkgTemplateMorphPdf_qqzz = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,SD,D,false,ROOT.RooArgList(bkgTemplatePdf_qqzz),ROOT.RooArgList(),1.0,1)
        TemplateName = "bkgTemplateMorphPdf_ggzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        bkgTemplateMorphPdf_ggzz = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,SD,D,false,ROOT.RooArgList(bkgTemplatePdf_ggzz),ROOT.RooArgList(),1.0,1)

        print 'MORPHING the ZJETS BKG'
        TemplateName = "bkgTemplateMorphPdf_zjets_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        bkgTemplateMorphPdf_zjets = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,SD,D,false,funcList_zjets,morphVarListBkg,1.0,1)

        print 'dummy cloning'
        bkg2d_qqzz =bkg_qqzz.Clone("bkg2d_qqzz")
        bkg2d_ggzz =bkg_ggzz.Clone("bkg2d_ggzz")
        bkg2d_zjets =bkg_zjets.Clone("bkg2d_zjets")
###        bkg2d_qqzz = ROOT.RooProdPdf("bkg2d_qqzz","bkg2d_qqzz",ROOT.RooArgSet(bkg_qqzz),ROOT.RooFit.Conditional(ROOT.RooArgSet(bkgTemplateMorphPdf_qqzz),ROOT.RooArgSet(D)))
###        bkg2d_ggzz = ROOT.RooProdPdf("bkg2d_ggzz","bkg2d_ggzz",ROOT.RooArgSet(bkg_ggzz),ROOT.RooFit.Conditional(ROOT.RooArgSet(bkgTemplateMorphPdf_ggzz),ROOT.RooArgSet(D)))
###        bkg2d_zjets = ROOT.RooProdPdf("bkg2d_zjets","bkg2d_zjets",ROOT.RooArgSet(bkg_zjets),ROOT.RooFit.Conditional(ROOT.RooArgSet(bkgTemplateMorphPdf_zjets),ROOT.RooArgSet(D)))

        ## ----------------- SUPERMELA 1D BACKGROUND SHAPES --------------- ##

        templateSDBkgName = "{0}/Dbackground_qqZZ_{1}.root".format(self.templateDir ,self.appendName) 
        print 'SUPERMELA 1D BACKGROUND SHAPES, ',templateSDBkgName
        bkgTempSDFile = ROOT.TFile(templateSDBkgName)
        bkgTemplateSD = bkgTempSDFile.Get("h_superD")
###        bkgTemplateSD = bkgTempSDFile.Get("h_superD_mod_rndm") 
###        bkgTemplateSD = bkgTempSDFile.Get("h_superDfromProjX") 
        TemplateSDName = "bkgTempSDDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        bkgTempSDDataHist = ROOT.RooDataHist(TemplateSDName,TemplateSDName,ROOT.RooArgList(SD),bkgTemplateSD)
        
        TemplateSDName = "bkgTemplateSDPdf_qqzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        bkgTemplateSDPdf_qqzz = ROOT.RooHistPdf(TemplateSDName,TemplateSDName,ROOT.RooArgSet(SD),bkgTempSDDataHist)

        TemplateSDName = "bkgTemplateSDPdf_ggzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        bkgTemplateSDPdf_ggzz = ROOT.RooHistPdf(TemplateSDName,TemplateSDName,ROOT.RooArgSet(SD),bkgTempSDDataHist)
        TemplateSDName = "bkgTemplateSDPdf_zjets_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        bkgTemplateSDPdf_zjets = ROOT.RooHistPdf(TemplateSDName,TemplateSDName,ROOT.RooArgSet(SD),bkgTempSDDataHist)
        

        ## ----------------------- PLOTS FOR SANITY CHECKS -------------------------- ##
     #    print 'PLOTS FOR SANITY CHECK'
#         czz = ROOT.TCanvas( "czz", "czz", 750, 700 )
#         czz.cd()
#         zzframe_s = CMS_zz4l_mass.frame(45)
#         super(RooDoubleCB,signalCB_ggH).plotOn(zzframe_s, ROOT.RooFit.LineStyle(3), ROOT.RooFit.LineColor(kYellow+3) )
#         super(ROOT.RooFFTConvPdf,sig_ggH).plotOn(zzframe_s, ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineColor(1) )
#         super(ROOT.RooqqZZPdf_v2,bkg_qqzz).plotOn(zzframe_s, ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineColor(4) )
#         super(ROOT.RooggZZPdf_v2,bkg_ggzz).plotOn(zzframe_s, ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineColor(6) )
#         super(ROOT.RooLandau,bkg_zjets).plotOn(zzframe_s, ROOT.RooFit.LineStyle(2), ROOT.RooFit.LineColor(6) )
#         zzframe_s.Draw()
#         figName = "{0}/figs/mzz_{1}_{2}.png".format(self.outputDir, self.mH, self.appendName)
#         czz.SaveAs(figName)
#         del czz
        
        ## ------------------- LUMI -------------------- ##
        print 'LUMI'
        rrvLumi = ROOT.RooRealVar("cmshzz4l_lumi","cmshzz4l_lumi",self.lumi)  
        
        ## ----------------------- SIGNAL RATES ----------------------- ##
        
        CMS_zz4l_mass.setRange("shape",self.low_M,self.high_M)
        
        fr_low_M = self.low_M
        fr_high_M = self.high_M        
        if (self.mH >= 450): 
            fr_low_M = 100
            fr_high_M = 1000
        if (self.mH >= 750):
            fr_low_M = 100
            fr_high_M = 1400
            
        print 'CMS_zz4l_mass range: ',fr_low_M,'  ',fr_high_M
        CMS_zz4l_mass.setRange("fullrangesignal",fr_low_M,fr_high_M)
        CMS_zz4l_mass.setRange("fullrange",100,1400)
        
        rfvCsFilter = RooFormulaVar()
        filterName = "cmshzz4l_csFilter_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.sqrts == 7): 
            rfvCsFilter = ROOT.RooFormulaVar(filterName,"0.5+0.5*TMath::Erf((@0 - 80.85)/50.42)", ROOT.RooArgList(self.MH) )
        else:
            rfvCsFilter = ROOT.RooFormulaVar(filterName,"@0",ROOT.RooArgList(one))

        # if(DEBUG):
        print "@@@@@@@ rfvCsFilter = ",rfvCsFilter.getVal()

        sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_a1".format(self.channel,self.sqrts)
        rrva1 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_a1'])
        sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_a2".format(self.channel,self.sqrts)
        rrva2 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_a2'])
        sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_a3".format(self.channel,self.sqrts)
        rrva3 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_a3'])
        sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_a4".format(self.channel,self.sqrts)
        rrva4 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_a4'])
        sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_b1".format(self.channel,self.sqrts)
        rrvb1 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_b1'])
        sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_b2".format(self.channel,self.sqrts)
        rrvb2 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_b2'])
        sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_b3".format(self.channel,self.sqrts)
        rrvb3 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_b3'])

        if(DEBUG):
            print "sigEff_a1 = ",theInputs['sigEff_a1']
            print "sigEff_a2 = ",theInputs['sigEff_a2']
            print "sigEff_a3 = ",theInputs['sigEff_a3']
            print "sigEff_a4 = ",theInputs['sigEff_a4']
            print "sigEff_b1 = ",theInputs['sigEff_b1']
            print "sigEff_b2 = ",theInputs['sigEff_b2']
            print "sigEff_b3 = ",theInputs['sigEff_b3']
           
    
        sigEffName = "hzz4lsignaleff_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)

        rfvSigEff = ROOT.RooFormulaVar(sigEffName,"(@0+@1*TMath::Erf((@7-@2)/@3))*(@4+@5*@7+@6*@7*@7)", ROOT.RooArgList(rrva1,rrva2,rrva3,rrva4,rrvb1,rrvb2,rrvb3,self.MH))
        ## following printout is needed ,  dont remove it
        print " @@@@@@@@ sigeff ",rfvSigEff.getVal()
    
        CS_ggH = myCSW.HiggsCS(1,self.mH,self.sqrts)
        CS_VBF = myCSW.HiggsCS(2,self.mH,self.sqrts)
        CS_WH = myCSW.HiggsCS(3,self.mH,self.sqrts)
        CS_ZH = myCSW.HiggsCS(4,self.mH,self.sqrts)
        CS_ttH = myCSW.HiggsCS(5,self.mH,self.sqrts)
    
        BRH2e2mu = myCSW.HiggsBR(13,self.mH)
        BRH4mu = myCSW.HiggsBR(12,self.mH)
        BRH4e = myCSW.HiggsBR(12,self.mH)
        BR = 0.0
        if( self.channel == self.ID_4mu ): BR = BRH4mu
        if( self.channel == self.ID_4e ): BR = BRH4e
        if( self.channel == self.ID_2e2mu ): BR = BRH2e2mu
    
        sigEfficiency = rfvSigEff.getVal()

        if(DEBUG):
            print "CS_ggH: ",CS_ggH,", CS_VBF: ",CS_VBF,", CS_WH: ",CS_WH,", CS_ZH: ",CS_ZH
            print ", CS_ttH: ",CS_ttH,", BRH2e2mu: ",BRH2e2mu,", BRH4mu: ",BRH4mu,", BRH4e: ",BRH4e

        csCorrection = 1.0
        if(self.sqrts == 7): csCorrection = self.csFilter(self.mH)

        ## SIG YIELDS
        sigRate_ggH = csCorrection*CS_ggH*BR*sigEfficiency*1000.*self.lumi
        sigRate_VBF = csCorrection*CS_VBF*BR*sigEfficiency*1000.*self.lumi
        sigRate_WH = csCorrection*CS_WH*BR*sigEfficiency*1000.*self.lumi
        sigRate_ZH = csCorrection*CS_ZH*BR*sigEfficiency*1000.*self.lumi
        sigRate_ttH = csCorrection*CS_ttH*BR*sigEfficiency*1000.*self.lumi
       
        normalizationSignal = signalCB_ggH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("fullrangesignal") ).getVal()
        
            
        print "#################### ",signalCB_ggH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("fullrangesignal") ).getVal()
        print "#################### ",signalCB_ggH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
   ####     print "#################### ",sig_ggH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("fullrangesignal") ).getVal()
        
        sclFactorSig_ggH = sigRate_ggH/normalizationSignal
        integral_ggH = 0.0
     
        integral_ggH =signalCB_ggH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        sigRate_ggH_Shape = sclFactorSig_ggH*integral_ggH
           
        normSigName = "cmshzz4l_normalizationSignal_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        rrvNormSig = ROOT.RooRealVar()


        rrvNormSig = ROOT.RooRealVar(normSigName,normSigName, self.getVariable(signalCB_ggH.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass)).getVal(),sig_ggH.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass)).getVal(),self.bUseCBnoConvolution))
        rrvNormSig.setConstant(True)

      
            
      #  rfvSigRate_ggH = ROOT.RooFormulaVar("ggH_norm","@0*@1*@2*1000*{0}*{2}/{1}".format(self.lumi,rrvNormSig.getVal(),self.getVariable(signalCB_ggH.createIntegral(RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal(),sig_ggH.createIntegral(RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal(),self.bUseCBnoConvolution)),ROOT.RooArgList(rfvCsFilter,rfvSigEff, rhfXsBrFuncV_1))

        rfvSigRate_ggH = ROOT.RooFormulaVar("ggH_norm","@0*@1*@2*1000*{0}*{2}/{1}".format(self.lumi,rrvNormSig.getVal(),integral_ggH),ROOT.RooArgList(rfvCsFilter,rfvSigEff, rhfXsBrFuncV_1))

        print "Compare integrals: integral_ggH=",integral_ggH,"  ; calculated=",self.getVariable(signalCB_ggH.createIntegral(RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal(),sig_ggH.createIntegral(RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal(),self.bUseCBnoConvolution)
        

      ###  rfvSigRate_ggH = ROOT.RooFormulaVar("ggH_norm","@0*@1*@2*1000*{0}*{2}/{1}".format(self.lumi,rrvNormSig.getVal(),integral_ggH),ROOT.RooArgList(rfvCsFilter,rfvSigEff, rhfXsBrFuncV_1))

        print "Compare integrals: integral_ggH=",integral_ggH,"  ; calculated=",self.getVariable(signalCB_ggH.createIntegral(RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal(),sig_ggH.createIntegral(RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal(),self.bUseCBnoConvolution)



        print signalCB_ggH.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass)).getVal(),"   ",sig_ggH.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass)).getVal()
        print signalCB_ggH.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal(),"   ",sig_ggH.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal()
        if (self.all_chan):
            print "Requested to sum up over the 5 chans: the norm in rfvSigRate_ggH should be the sum of the values of sigRate_XYZ_Shape variables:"
        print " @@@@@@@ norm sig = ",rrvNormSig.getVal()
        print " @@@@@@@ rfvSigRate_ggH = ",rfvSigRate_ggH.getVal()
        print " sigRate_ggH_Shape=",sigRate_ggH_Shape
   
   #     print "Sum of sigRate_XYZ_Shape=",sigRate_ggH_Shape+sigRate_VBF_Shape+sigRate_WH_Shape+sigRate_ZH_Shape+sigRate_ttH_Shape
        ## SET RATES TO 1 
        ## DC RATES WILL BE MULTIPLIED
        ## BY RATES IMPORTED TO WS
        sigRate_ggH_Shape = 1
        sigRate_VBF_Shape = 1
        sigRate_WH_Shape = 1
        sigRate_ZH_Shape = 1
        sigRate_ttH_Shape = 1

             
        ## ----------------------- BACKGROUND RATES ----------------------- ##

        ## rates per lumi for scaling
        bkgRate_qqzz = theInputs['qqZZ_rate']/theInputs['qqZZ_lumi']
        bkgRate_ggzz = theInputs['ggZZ_rate']/theInputs['qqZZ_lumi']
        bkgRate_zjets = theInputs['zjets_rate']/theInputs['zjets_lumi']
        
        ## Get Normalizations
        normalizationBackground_qqzz = bkg_qqzz.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("fullrange") ).getVal()
        normalizationBackground_ggzz = bkg_ggzz.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("fullrange") ).getVal()
        normalizationBackground_zjets = bkg_zjets.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("fullrange") ).getVal()
        
        sclFactorBkg_qqzz = self.lumi*bkgRate_qqzz/normalizationBackground_qqzz
        sclFactorBkg_ggzz = self.lumi*bkgRate_ggzz/normalizationBackground_ggzz
        sclFactorBkg_zjets = self.lumi*bkgRate_zjets/normalizationBackground_zjets
               
        bkgRate_qqzz_Shape = sclFactorBkg_qqzz * bkg_qqzz.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        bkgRate_ggzz_Shape = sclFactorBkg_ggzz * bkg_ggzz.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        bkgRate_zjets_Shape = sclFactorBkg_zjets * bkg_zjets.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        
        if(DEBUG):
            print "Shape signal rate: ",sigRate_ggH_Shape,", background rate: ",bkgRate_qqzz_Shape,", ",bkgRate_zjets_Shape," in ",self.low_M," - ",self.high_M
            CMS_zz4l_mass.setRange("lowmassregion",100.,160.)
            bkgRate_qqzz_lowmassregion = sclFactorBkg_qqzz * bkg_qqzz.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("lowmassregion") ).getVal()
            bkgRate_ggzz_lowmassregion = sclFactorBkg_ggzz * bkg_ggzz.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("lowmassregion") ).getVal()
            bkgRate_zjets_lowmassregion = sclFactorBkg_zjets * bkg_zjets.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("lowmassregion") ).getVal()
            lowmassyield = bkgRate_qqzz_lowmassregion + bkgRate_ggzz_lowmassregion + bkgRate_zjets_lowmassregion
            print "low mass yield: ",lowmassyield
        
        ## --------------------------- DATASET --------------------------- ##

        dataFileDir = "CMSdata"
        dataTreeName = "data_obs"  ### "data_obs" "SelectedTree" 
        dataFileName = "{0}/hzz{1}_{2}.root".format(dataFileDir,self.appendName,self.lumi)
###        dataFileName = "{0}/hzzDummyTMP_{1}_{2}_withSMD_doubleCBonly.root".format(dataFileDir,self.appendName,self.lumi)
        if (DEBUG): print dataFileName," ",dataTreeName 
        data_obs_file = ROOT.TFile(dataFileName)

        print data_obs_file.Get(dataTreeName)
        
        if not (data_obs_file.Get(dataTreeName)):
            print "File, \"",dataFileName,"\", or tree, \"",dataTreeName,"\", not found" 
            print "Exiting..."
            sys.exit()
        
        data_obs_tree = data_obs_file.Get(dataTreeName)
        data_obs = ROOT.RooDataSet()
        datasetName = "data_obs_{0}".format(self.appendName)
        
###            data_obs = ROOT.RooDataSet(datasetName,datasetName,data_obs_tree,ROOT.RooArgSet(SD,D))
        data_obs = ROOT.RooDataSet(datasetName,datasetName,data_obs_tree,ROOT.RooArgSet(CMS_zz4l_mass,SD,D),'CMS_zz4l_mass>106.0&&CMS_zz4l_mass<141.0').reduce(ROOT.RooArgSet(SD,D))


            
        ## --------------------------- WORKSPACE -------------------------- ##

        endsInP5 = False
        tmpMH = self.mH
        if ( math.fabs(math.floor(tmpMH)-self.mH) == 0.5): endsInP5 = True
        if (DEBUG): print "ENDS IN P5  ",endsInP5

        name_Shape = ""
        name_ShapeWS = ""
        name_ShapeWS2 = ""
        name_ShapeWSXSBR = ""
        
        if (endsInP5): name_Shape = "{0}/HCG/{1}/hzz4l_{2}S_{3:.0f}TeV.txt".format(self.outputDir,self.mH,self.appendName,self.sqrts)
        else: name_Shape = "{0}/HCG/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV.txt".format(self.outputDir,self.mH,self.appendName,self.sqrts)

        if (endsInP5): name_ShapeWS = "{0}/HCG/{1}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts)
        else: name_ShapeWS = "{0}/HCG/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts)

        if (endsInP5): name_ShapeWSXSBR = "{0}/HCG_XSxBR/{1}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts)
        else: name_ShapeWSXSBR = "{0}/HCG_XSxBR/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts)

        name_ShapeWS2 = "hzz4l_{0}S_{1:.0f}TeV.input.root".format(self.appendName,self.sqrts)

        if(DEBUG): print name_Shape,"  ",name_ShapeWS2
        
        w = ROOT.RooWorkspace("w","w")
        
        w.importClassCode(RooqqZZPdf_v2.Class(),True)
        w.importClassCode(RooggZZPdf_v2.Class(),True)
        w.importClassCode(RooRelBWUFParam.Class(),True)
        w.importClassCode(RooDoubleCB.Class(),True)
        w.importClassCode(RooFormulaVar.Class(),True)


        if( FactorizedShapes ):
            if( self.channel == self.ID_4mu ):
                w.importClassCode(RooFourMuMassShapePdf2.Class(),True)
                w.importClassCode(RooFourMuMassRes.Class(),True)
            elif( self.channel == self.ID_4e ):
                w.importClassCode(RooFourEMassShapePdf2.Class(),True)
                w.importClassCode(RooFourEMassRes.Class(),True)
            elif( self.channel == self.ID_2e2mu ):
                w.importClassCode(RooTwoETwoMuMassShapePdf2.Class(),True)
                w.importClassCode(RooTwoETwoMuMassRes.Class(),True)
            
                
                
        getattr(w,'import')(data_obs,ROOT.RooFit.Rename("data_obs")) ### Should this be renamed?
    
        if(self.bUseCBnoConvolution) :
            print 'Saving objects in WS (bUseCBnoConvolution==true)'

###################################
#### WE USE THIS SNIPPET ##########                                         

            if (self.is2D == 1):
                print 'Saving 2D templates with SuperMELA in the WS (bUseCBnoConvolution==true)'
                sigCB2d_ggH.SetNameTitle("ORIGggH","ORIGggH")
               
                getattr(w,'import')(sigCB2d_ggH, ROOT.RooFit.RecycleConflictNodes())

                #  sigTemplatePdf_ggH.SetNameTitle("ggH","ggH")
                sigTemplateMorphPdf_ggH.SetNameTitle("ggH","ggH")
                getattr(w,'import')(sigTemplateMorphPdf_ggH, ROOT.RooFit.RecycleConflictNodes())
                #save syst templates individually
                systTempName=("ggHCMS_zz4l_leptScale_sig_{0}_{1:.0f}_Up").format(self.channel,self.sqrts)
                sigTemplatePdf_ggH_syst1Up.SetNameTitle(systTempName,systTempName)
                systTempName=("ggHCMS_zz4l_leptScale_sig_{0}_{1:.0f}_Down").format(self.channel,self.sqrts)
                sigTemplatePdf_ggH_syst1Down.SetNameTitle(systTempName,systTempName)
                
                if(self.isAltSig):   
                    sigCB2d_ggH_ALT.SetNameTitle("ORIGggH{0}".format(self.appendHypType),"ggH{0}".format(self.appendHypType))
                    #       sigTemplatePdf_ggH_ALT.SetNameTitle("ggH{0}".format(self.appendHypType),"ggH{0}".format(self.appendHypType))
                    sigTemplateMorphPdf_ggH_ALT.SetNameTitle("ggH{0}".format(self.appendHypType),"ggH{0}".format(self.appendHypType))
                    getattr(w,'import')(sigCB2d_ggH_ALT, ROOT.RooFit.RecycleConflictNodes())
                    
                    getattr(w,'import')(sigTemplatePdf_ggH_ALT, ROOT.RooFit.RecycleConflictNodes())
                    getattr(w,'import')(sigTemplateMorphPdf_ggH_ALT, ROOT.RooFit.RecycleConflictNodes())
                    

            print 'Saving in WS bkgd shapes for superMELA'
            bkg2d_qqzz.SetNameTitle("ORIGbkg2d_qqzz","ORIGbkg2d_qqzz")
            bkg2d_ggzz.SetNameTitle("ORIGbkg2d_ggzz","ORIGbkg2d_ggzz")
            bkg2d_zjets.SetNameTitle("ORIGbkg2d_zjets","ORIGbkg2d_zjets")            
    
            bkgTemplateMorphPdf_qqzz.SetNameTitle("bkg2d_qqzz","bkg2d_qqzz")
            bkgTemplateMorphPdf_ggzz.SetNameTitle("bkg2d_ggzz","bkg2d_ggzz")

            bkgTemplateMorphPdf_zjets.SetNameTitle("bkg2d_zjets","bkg2d_zjets")
            
            getattr(w,'import')(bkgTemplateMorphPdf_qqzz,ROOT.RooFit.RecycleConflictNodes())
            getattr(w,'import')(bkgTemplateMorphPdf_ggzz,ROOT.RooFit.RecycleConflictNodes())
            getattr(w,'import')(bkgTemplateMorphPdf_zjets,ROOT.RooFit.RecycleConflictNodes())
            print 'finished to save bkg 2D pdfs in the WS'
#################


### save signal rates        
        getattr(w,'import')(rfvSigRate_ggH, ROOT.RooFit.RecycleConflictNodes())
        
        if(self.isAltSig):

            rfvSigRate_ggH_ALT=ROOT.RooFormulaVar(rfvSigRate_ggH,"ggH{0}_norm".format(self.appendHypType))
            print 'Compare signal rates: STD=',rfvSigRate_ggH.getVal(),"   ALT=",rfvSigRate_ggH_ALT.getVal()
            getattr(w,'import')(rfvSigRate_ggH_ALT, ROOT.RooFit.RecycleConflictNodes())
            
        w.writeToFile(name_ShapeWS)
        w.writeToFile(name_ShapeWSXSBR)
        
        ## --------------------------- DATACARDS -------------------------- ##

        systematics.setSystematics(bkgRate_qqzz_Shape, bkgRate_ggzz_Shape, bkgRate_zjets_Shape)
        systematics_forXSxBR.setSystematics(bkgRate_qqzz_Shape, bkgRate_ggzz_Shape,bkgRate_zjets_Shape)

        ## If the channel is not declared in inputs, set rate = 0
        if not self.ggH_chan and not self.all_chan :  sigRate_ggH_Shape = 0
        if not self.qqH_chan:  sigRate_VBF_Shape = 0
        if not self.WH_chan:   sigRate_WH_Shape = 0
        if not self.ZH_chan:   sigRate_ZH_Shape = 0
        if not self.ttH_chan:  sigRate_ttH_Shape = 0

        if not self.qqZZ_chan:  bkgRate_qqzz_Shape = 0
        if not self.ggZZ_chan:  bkgRate_ggzz_Shape = 0
        if not self.zjets_chan: bkgRate_zjets_Shape = 0

        rates = {}
        rates['ggH'] = sigRate_ggH_Shape
        rates['qqH'] = sigRate_VBF_Shape
        rates['WH']  = sigRate_WH_Shape
        rates['ZH']  = sigRate_ZH_Shape
        rates['ttH'] = sigRate_ttH_Shape

        rates['qqZZ']  = bkgRate_qqzz_Shape
        rates['ggZZ']  = bkgRate_ggzz_Shape
        rates['zjets'] = bkgRate_zjets_Shape
        rates['ttbar'] = 0
        rates['zbb']   = 0
        

        ## Write Datacards
        fo = open( name_Shape, "wb")
        self.WriteDatacard(fo,theInputs, name_ShapeWS2, rates, data_obs.numEntries(), self.is2D )
        systematics.WriteSystematics(fo, theInputs,)
        systematics.WriteShapeSystematics(fo,theInputs)
        fo.close()

        if(self.isAltSig):
            
            if (endsInP5): name_Shape = "{0}/HCG/{1}/hzz4l_{2}S_{3:.0f}TeV{4}.txt".format(self.outputDir,self.mH,self.appendName,self.sqrts,self.appendHypType)
            else: name_Shape = "{0}/HCG/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV{4}.txt".format(self.outputDir,self.mH,self.appendName,self.sqrts,self.appendHypType)
            fo = open( name_Shape, "wb")
            self.WriteDatacard(fo,theInputs, name_ShapeWS2, rates, data_obs.numEntries(), self.is2D,True,self.appendHypType )
            systematics.WriteSystematics(fo, theInputs,self.isAltSig)
            systematics.WriteShapeSystematics(fo,theInputs)
            fo.close()

        ## forXSxBR
        if (endsInP5): name_Shape = "{0}/HCG_XSxBR/{2}/hzz4l_{1}S_{3:.0f}TeV.txt".format(self.outputDir,self.appendName,self.mH,self.sqrts)	
        else: name_Shape = "{0}/HCG_XSxBR/{2:.0f}/hzz4l_{1}S_{3:.0f}TeV.txt".format(self.outputDir,self.appendName,self.mH,self.sqrts)
            
        fo = open( name_Shape, "wb" )
        self.WriteDatacard(fo, theInputs,name_ShapeWS2, rates, data_obs.numEntries(), self.is2D )
        systematics_forXSxBR.WriteSystematics(fo, theInputs)
        systematics_forXSxBR.WriteShapeSystematics(fo,theInputs)
        fo.close()

        if(self.isAltSig):
            if (endsInP5): name_Shape = "{0}/HCG_XSxBR/{2}/hzz4l_{1}S_{3:.0f}TeV{4}.txt".format(self.outputDir,self.appendName,self.mH,self.sqrts,self.appendHypType)	
            else: name_Shape = "{0}/HCG_XSxBR/{2:.0f}/hzz4l_{1}S_{3:.0f}TeV{4}.txt".format(self.outputDir,self.appendName,self.mH,self.sqrts,self.appendHypType)
            fo = open( name_Shape, "wb")
            self.WriteDatacard(fo,theInputs,name_ShapeWS2,rates,data_obs.numEntries(),self.is2D,True,self.appendHypType )
            systematics.WriteSystematics(fo, theInputs,self.isAltSig)
            systematics.WriteShapeSystematics(fo,theInputs)
            fo.close()
        



    def WriteDatacard(self,file,theInputs,nameWS,theRates,obsEvents,is2D,isAltCard=False,AltLabel=""):

        numberSig = self.numberOfSigChan(theInputs)
        numberBg  = self.numberOfBgChan(theInputs)
        if isAltCard and theInputs["all"] : numberSig = 2
        
        file.write("imax 1\n")
        file.write("jmax {0}\n".format(numberSig+numberBg-1))
        file.write("kmax *\n")
        
        file.write("------------\n")
        file.write("shapes * * {0} w:$PROCESS \n".format(nameWS))
        file.write("------------\n")
        
        file.write("bin a{0} \n".format(self.channel))
        file.write("observation {0} \n".format(obsEvents))
        
        file.write("------------\n")
        file.write("## mass window [{0},{1}] \n".format(self.low_M,self.high_M))
        file.write("bin ")        

        channelList=['ggH','qqH','WH','ZH','ttH','qqZZ','ggZZ','zjets','ttbar','zbb']

        channelName1D=['ggH','qqH','WH','ZH','ttH','bkg_qqzz','bkg_ggzz','bkg_zjets','bkg_ttbar','bkg_zbb']
        channelName2D=['ggH','qqH','WH','ZH','ttH','bkg2d_qqzz','bkg2d_ggzz','bkg2d_zjets','bkg2d_ttbar','bkg2d_zbb']

        if theInputs["all"]:
           channelList=['ggH','qqZZ','ggZZ','zjets','ttbar','zbb']
           channelName1D=['ggH','bkg_qqzz','bkg_ggzz','bkg_zjets','bkg_ttbar','bkg_zbb']
           channelName2D=['ggH','bkg2d_qqzz','bkg2d_ggzz','bkg2d_zjets','bkg2d_ttbar','bkg2d_zbb']

           if isAltCard :
               print 'THIS IS AN ALTERNATIVE CARD !!!!'
               channelList=['ggH','ggH','qqZZ','ggZZ','zjets','ttbar','zbb']
               channelName1D=['ggH','ggH{0}'.format(AltLabel),'bkg_qqzz','bkg_ggzz','bkg_zjets','bkg_ttbar','bkg_zbb']
               channelName2D=['ggH','ggH{0}'.format(AltLabel),'bkg2d_qqzz','bkg2d_ggzz','bkg2d_zjets','bkg2d_ttbar','bkg2d_zbb']



        if ( len(channelList) != len(channelName2D) ) :
            raise RuntimeError, "Mismatch in length of channel arrays!"

        for chan in channelList:
            if theInputs[chan]:
                file.write("a{0} ".format(self.channel))
            else:
                if chan.startswith("ggH") and theInputs["all"] :
                    file.write("a{0} ".format(self.channel))
        file.write("\n")
                                        
        file.write("process ")

        i=0

        for chan in channelList:
            #  print 'checking if ',chan,' is in the list of to-do'
            if theInputs[chan]:
                file.write("{0} ".format(channelName2D[i]))
                #  print 'writing in card index=',i,'  chan=',chan
                i+=1
            else:
                if chan.startswith("ggH") and theInputs["all"] :
                    file.write("{0} ".format(channelName2D[i]))
                    # print 'writing in card TOTAL SUM, index=',i,'  chan=',chan,'  ',channelName2D[i]
                    i+=1
        
        file.write("\n")
            
        processLine = "process "

        for x in range(-numberSig+1,1):
            processLine += "{0} ".format(x)

        for y in range(1,numberBg+1):
            processLine += "{0} ".format(y)

        file.write(processLine)
        file.write("\n")
            
        file.write("rate ")
        ### correction to the yield related to interference effects
        ### (depends on spin hypothesis, only for 4e and 4mu)
   ##     corrForInterf = self.calcInterfYieldCorr(self.channel,self.spinHyp)
   ##     corrForAccept = self.calcAcceptYieldCorr(self.channel,self.spinHyp)
        corrForInterfAndAccept = self.calcTotalYieldCorr(self.channel,self.spinHyp)
   ###     print 'Corr for interference is ',corrForInterf
        ind=0
        for chan in channelList:
            if theInputs[chan] or (chan.startswith("ggH") and theInputs["all"]):
      ###          print 'ChannelName2D=',channelName2D[ind]
                if ( channelName2D[ind]=="ggH_ALT") :
   ###                 print 'I will write in the card a factor ',theRates[chan]*corrForInterf
###                    file.write("{0:.4f} ".format(theRates[chan]*corrForInterf*corrForAccept))
                    file.write("{0:.4f} ".format(theRates[chan]*corrForInterfAndAccept))
                else :
                    file.write("{0:.4f} ".format(theRates[chan]))
            ind += 1
        file.write("\n")
        file.write("------------\n")

        
    def numberOfSigChan(self,inputs):

        counter=0

        if inputs['ggH']: counter+=1
        if inputs['qqH']: counter+=1
        if inputs['WH']:  counter+=1
        if inputs['ZH']:  counter+=1
        if inputs['ttH']: counter+=1
        if inputs['all']: counter+=1
        
        return counter

    def numberOfBgChan(self,inputs):

        counter=0

        if inputs['qqZZ']:  counter+=1
        if inputs['ggZZ']:  counter+=1
        if inputs['zjets']: counter+=1
        if inputs['ttbar']: counter+=1
        if inputs['zbb']:   counter+=1
        
        return counter

    def calcInterfYieldCorr(self, channel, spinHypCode):

        ### the parameters are the fraction of 2e2mu with respect to the total
        ### we modify the total assuming that the 2e2mu yield is constant
        ### among different spin hypotheses

        cSM=0.4816
        cALT=0.4816 ### do nothing by default

        if (spinHypCode == 0) : ### 0+ vs 0-
            cALT = 0.5236
        elif (spinHypCode == 1 or spinHypCode == 5) : ### 0+ vs 2m+ from gg or qq
            cALT = 0.5265 #does not depend on production mechanism
        elif (spinHypCode == 2) : ### 0+ vs 0h+
            cALT = 0.5084
        elif (spinHypCode == 3) : ### 0+ vs 1+
            cALT = 0.5210
        elif (spinHypCode == 4) : ### 0+ vs 1-
            cALT = 0.5068
        elif (spinHypCode < 0) : ### dummy/default
            cALT = cSM
        else :
            print 'Unrecognized SpinHypothesis=',spinHypCode,'  will not correct total 4e and 4mu yields'

        corrFactor=0.0
        if( channel == self.ID_2e2mu ):
            corrFactor=1.0
        else :
            corrFactor= ( cSM*(1.0-cALT) ) / (cALT*(1.0-cSM))
            
        print 'calcInterfYieldCorr: HypCode=',spinHypCode,'  cALT=',cALT,'  corrFactor=',corrFactor
        return corrFactor


            

    def calcAcceptYieldCorr(self, channel, spinHypCode):

        ### the parameters are the fraction of 2e2mu with respect to the total
        ### we modify the total assuming that the 2e2mu yield is constant
        ### among different spin hypotheses


        alpha=1.0 ### ratio of generated events w.r.t. SM 0+ sample
        r=1.0 ### ratio of reconstructed events w.r.t. SM0+

        if (spinHypCode == 0) : ### 0+ vs 0-
            if( channel == self.ID_2e2mu ):
                alpha=1.0
                r=1.0
            elif ( channel == self.ID_4mu ):
                alpha=1.0
                r=1.0
            else :  #( channel == self.ID_4e )
                alpha=1.0
                r=1.0
        elif (spinHypCode == 1) : ### gg->0+ vs 2m+
            if( channel == self.ID_2e2mu ):
                alpha=1.0
                r=1.0
            elif ( channel == self.ID_4mu ):
                alpha=1.0
                r=1.0
            else :  #( channel == self.ID_4e )
                alpha=1.0
                r=1.0
        elif (spinHypCode == 2) : ###0+ vs 0h+
            if( channel == self.ID_2e2mu ):
                alpha=1.0
                r=1.0
            elif ( channel == self.ID_4mu ):
                alpha=1.0
                r=1.0
            else :  #( channel == self.ID_4e )
                alpha=1.0
                r=1.0
        elif (spinHypCode == 3) : ### 0+ vs 1+
            if( channel == self.ID_2e2mu ):
                alpha=1.0
                r=1.0
            elif ( channel == self.ID_4mu ):
                alpha=1.0
                r=1.0
            else :  #( channel == self.ID_4e )
                alpha=1.0
                r=1.0
        elif (spinHypCode == 4) : ### 0+ vs 1-
            if( channel == self.ID_2e2mu ):
                alpha=1.0
                r=1.0
            elif ( channel == self.ID_4mu ):
                alpha=1.0
                r=1.0
            else :  #( channel == self.ID_4e )
                alpha=1.0
                r=1.0
        elif (spinHypCode == 5) : ### qq->0+ vs 2m+
            if( channel == self.ID_2e2mu ):
                alpha=1.0
                r=1.0
            elif ( channel == self.ID_4mu ):
                alpha=1.0
                r=1.0
            else :  #( channel == self.ID_4e )
                alpha=1.0
                r=1.0
        elif (spinHypCode < 0) : ### do nothing
            if( channel == self.ID_2e2mu ):
                alpha=1.0
                r=1.0
            elif ( channel == self.ID_4mu ):
                alpha=1.0
                r=1.0
            else :  #( channel == self.ID_4e )
                alpha=1.0
                r=1.0
        else :
            print 'Unrecognized SpinHypothesis=',spinHypCode,'  will not correct total yields because of acceptance'

        corrFactor=0.0
        if( channel == self.ID_2e2mu ):
            corrFactor=1.0
        else :
            corrFactor= r / alpha
            
        print 'calcAcceptYieldCorr: HypCode=',spinHypCode,'  alpha=',alpha,' r=',r,'  corrFactor=',corrFactor
        return corrFactor



    def calcTotalYieldCorr(self, channel, spinHypCode):

        ### the parameters are the fraction of 2e2mu with respect to the total
        ### we modify the total assuming that the 2e2mu yield is constant
        ### among different spin hypotheses



        r=1.0 ### correction calculated by Andrew, embeds both Interf and Accept

        if (spinHypCode == 0) : ### 0+ vs 0-
            if( channel == self.ID_2e2mu ):
                r=1.12188
            elif ( channel == self.ID_4mu ):
                r= 0.912885
            else :  #( channel == self.ID_4e )
                r= 0.858691
        elif (spinHypCode == 1) : ### gg->0+ vs 2m+
            if( channel == self.ID_2e2mu ):
                r= 1.12826
            elif ( channel == self.ID_4mu ):
                r=0.900717
            else :  #( channel == self.ID_4e )
                r= 0.865112
        elif (spinHypCode == 2) : ###0+ vs 0h+
            if( channel == self.ID_2e2mu ):
                r=1.06477
            elif ( channel == self.ID_4mu ):
                r= 0.941477
            else :  #( channel == self.ID_4e )
                r= 0.947087
        elif (spinHypCode == 3) : ### 0+ vs 1+
            if( channel == self.ID_2e2mu ):
                r= 1.07114
            elif ( channel == self.ID_4mu ):
                r=0.968514
            else :  #( channel == self.ID_4e )
                r= 0.882397
        elif (spinHypCode == 4) : ### 0+ vs 1-
            if( channel == self.ID_2e2mu ):
                r=1.1036
            elif ( channel == self.ID_4mu ):
                r=0.939779
            else :  #( channel == self.ID_4e )
                r=0.854791
        elif (spinHypCode == 5) : ### qq->0+ vs 2m+
            if( channel == self.ID_2e2mu ):
                r=1.1417
            elif ( channel == self.ID_4mu ):
                r= 0.884544
            else :  #( channel == self.ID_4e )
                r= 0.861437
        elif (spinHypCode < 0) : ### do nothing
            if( channel == self.ID_2e2mu ):
                alpha=1.0
                r=1.0
            elif ( channel == self.ID_4mu ):
                alpha=1.0
                r=1.0
            else :  #( channel == self.ID_4e )
                alpha=1.0
                r=1.0
        else :
            print 'Unrecognized SpinHypothesis=',spinHypCode,'  will not correct total yields because of acceptance'

            
        print 'calcTotalYieldCorr: ',r
        return r 


            
