#! /usr/bin/env python
import sys, os, pwd, commands
import optparse, shlex, re
import math
from scipy.special import erf
from ROOT import *
import ROOT
from array import array
from systematicsClassSMD import *
from inputReader import *

##-------------------------------
## script to turn 2D templates into 1D
## script also makes datatrees compatable with these templates
##
##
##
##
## Written by Chris Martin
## based on HZZ4L_Combination/CombinationPy code

#function for parsing options
def parseOptions():

    usage = ('usage: %prog [options] datasetList\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    parser.add_option('-t', '--templateDir', type='string', dest='templateDir', default="templates2D" ,help='directory with 2D template histos')
    parser.add_option('-o', '--OutDir', type='string', dest='OutDir', default="templates1D", help='directory with 1D template histos')
    parser.add_option('-l', '--Lumi', type='string', dest='lumi',default="5.015", help='lumi')
    
    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()


# the main procedure
def make1Dtemplates():
    
    # parse the arguments and options
    global opt, args
    parseOptions()
    c=0
    while (c < 3 ):
        appendName = ''
        if (c == 0): appendName = '4mu'
        elif (c == 1): appendName = '4e'
        elif (c == 2): appendName = '2e2mu'

    
        templateSigName = "{0}/Dsignal_{1}.root".format(opt.templateDir,appendName)
        sigTempFile = ROOT.TFile(templateSigName)
        sigTemplate = sigTempFile.Get("h_superDpsD")
        sigTemplate_syst1Up = sigTempFile.Get("h_superDpsD_LeptScaleUp")
        sigTemplate_syst1Down = sigTempFile.Get("h_superDpsD_LeptScaleDown")
        sigTemplate_syst2Up = sigTempFile.Get("h_superDpsD_LeptSmearUp")
        sigTemplate_syst2Down = sigTempFile.Get("h_superDpsD_LeptSmearDown")

        outTemplate = ROOT.TH1F("h_superDpsD","h_superDpsD",(sigTemplate.GetNbinsX())*(sigTemplate.GetNbinsY()),0,(sigTemplate.GetNbinsX())*(sigTemplate.GetNbinsY()))
        outTemplate_syst1Up = ROOT.TH1F("h_superDpsD_LeptScaleUp","h_superDpsD_LeptScaleUp",(sigTemplate_syst1Up.GetNbinsX())*(sigTemplate_syst1Up.GetNbinsY()),0,(sigTemplate_syst1Up.GetNbinsX())*(sigTemplate_syst1Up.GetNbinsY()))
        outTemplate_syst1Down = ROOT.TH1F("h_superDpsD_LeptScaleDown","h_superDpsD_LeptScaleDown",(sigTemplate_syst1Down.GetNbinsX())*(sigTemplate_syst1Down.GetNbinsY()),0,(sigTemplate_syst1Down.GetNbinsX())*(sigTemplate_syst1Down.GetNbinsY()))
        outTemplate_syst2Up = ROOT.TH1F("h_superDpsD_LeptSmearUp","h_superDpsD_LeptSmearUp",(sigTemplate_syst2Up.GetNbinsX())*(sigTemplate_syst2Up.GetNbinsY()),0,(sigTemplate_syst2Up.GetNbinsX())*(sigTemplate_syst2Up.GetNbinsY()))
        outTemplate_syst2Down = ROOT.TH1F("h_superDpsD_LeptSmearDown","h_superDpsD_LeptSmearDown",(sigTemplate_syst2Down.GetNbinsX())*(sigTemplate_syst2Down.GetNbinsY()),0,(sigTemplate_syst2Down.GetNbinsX())*(sigTemplate_syst2Down.GetNbinsY()))

        bin_tot=0
        for xbin in xrange(sigTemplate.GetNbinsX()):
            for ybin in xrange(sigTemplate.GetNbinsY()):
                outTemplate.SetBinContent(bin_tot+1,sigTemplate.GetBinContent(xbin+1,ybin+1))
                outTemplate_syst1Up.SetBinContent(bin_tot+1,sigTemplate_syst1Up.GetBinContent(xbin+1,ybin+1))
                outTemplate_syst1Down.SetBinContent(bin_tot+1,sigTemplate_syst1Down.GetBinContent(xbin+1,ybin+1))
                outTemplate_syst2Up.SetBinContent(bin_tot+1,sigTemplate_syst2Up.GetBinContent(xbin+1,ybin+1))
                outTemplate_syst2Down.SetBinContent(bin_tot+1,sigTemplate_syst2Down.GetBinContent(xbin+1,ybin+1))
                bin_tot+=1

        OutTempSigName = "{0}/Dsignal_{1}.root".format(opt.OutDir,appendName)
        OutTempFile = ROOT.TFile(OutTempSigName,"RECREATE")
        outTemplate.Write()
        outTemplate_syst1Up.Write()
        outTemplate_syst1Down.Write()
        outTemplate_syst2Up.Write()
        outTemplate_syst2Down.Write()
        OutTempFile.Write()

        sigTempFile.Close()
        OutTempFile.Close()

        templateSigName_ALT = "{0}/Dsignal_ALT_{1}.root".format(opt.templateDir,appendName)
        sigTempFile_ALT = ROOT.TFile(templateSigName_ALT)
        sigTemplate_ALT = sigTempFile_ALT.Get("h_superDpsD")
        sigTemplate_ALT_syst1Up = sigTempFile_ALT.Get("h_superDpsD_LeptScaleUp")
        sigTemplate_ALT_syst1Down = sigTempFile_ALT.Get("h_superDpsD_LeptScaleDown")
        sigTemplate_ALT_syst2Up = sigTempFile_ALT.Get("h_superDpsD_LeptSmearUp")
        sigTemplate_ALT_syst2Down = sigTempFile_ALT.Get("h_superDpsD_LeptSmearDown")

        outTemplate_ALT = ROOT.TH1F("h_superDpsD","h_superDpsD",(sigTemplate_ALT.GetNbinsX())*(sigTemplate_ALT.GetNbinsY()),0,(sigTemplate_ALT.GetNbinsX())*(sigTemplate_ALT.GetNbinsY()))
        outTemplate_ALT_syst1Up = ROOT.TH1F("h_superDpsD_LeptScaleUp","h_superDpsD_LeptScaleUp",(sigTemplate_ALT_syst1Up.GetNbinsX())*(sigTemplate_ALT_syst1Up.GetNbinsY()),0,(sigTemplate_ALT_syst1Up.GetNbinsX())*(sigTemplate_ALT_syst1Up.GetNbinsY()))
        outTemplate_ALT_syst1Down = ROOT.TH1F("h_superDpsD_LeptScaleDown","h_superDpsD_LeptScaleDown",(sigTemplate_ALT_syst1Down.GetNbinsX())*(sigTemplate_ALT_syst1Down.GetNbinsY()),0,(sigTemplate_ALT_syst1Down.GetNbinsX())*(sigTemplate_ALT_syst1Down.GetNbinsY()))
        outTemplate_ALT_syst2Up = ROOT.TH1F("h_superDpsD_LeptSmearUp","h_superDpsD_LeptSmearUp",(sigTemplate_ALT_syst2Up.GetNbinsX())*(sigTemplate_ALT_syst2Up.GetNbinsY()),0,(sigTemplate_ALT_syst2Up.GetNbinsX())*(sigTemplate_ALT_syst2Up.GetNbinsY()))
        outTemplate_ALT_syst2Down = ROOT.TH1F("h_superDpsD_LeptSmearDown","h_superDpsD_LeptSmearDown",(sigTemplate_ALT_syst2Down.GetNbinsX())*(sigTemplate_ALT_syst2Down.GetNbinsY()),0,(sigTemplate_ALT_syst2Down.GetNbinsX())*(sigTemplate_ALT_syst2Down.GetNbinsY()))

        bin_ALT_tot=0
        for xbin in xrange(sigTemplate_ALT.GetNbinsX()):
            for ybin in xrange(sigTemplate_ALT.GetNbinsY()):
                outTemplate_ALT.SetBinContent(bin_ALT_tot+1,sigTemplate_ALT.GetBinContent(xbin+1,ybin+1))
                outTemplate_ALT_syst1Up.SetBinContent(bin_ALT_tot+1,sigTemplate_ALT_syst1Up.GetBinContent(xbin+1,ybin+1))
                outTemplate_ALT_syst1Down.SetBinContent(bin_ALT_tot+1,sigTemplate_ALT_syst1Down.GetBinContent(xbin+1,ybin+1))
                outTemplate_ALT_syst2Up.SetBinContent(bin_ALT_tot+1,sigTemplate_ALT_syst2Up.GetBinContent(xbin+1,ybin+1))
                outTemplate_ALT_syst2Down.SetBinContent(bin_ALT_tot+1,sigTemplate_ALT_syst2Down.GetBinContent(xbin+1,ybin+1))
                bin_ALT_tot+=1

        OutTempSigName_ALT = "{0}/Dsignal_ALT_{1}.root".format(opt.OutDir,appendName)
        OutTempFile_ALT = ROOT.TFile(OutTempSigName_ALT,"RECREATE")
        outTemplate_ALT.Write()
        outTemplate_ALT_syst1Up.Write()
        outTemplate_ALT_syst1Down.Write()
        outTemplate_ALT_syst2Up.Write()
        outTemplate_ALT_syst2Down.Write()
        OutTempFile_ALT.Write()

        sigTempFile_ALT.Close()
        OutTempFile_ALT.Close()

        templateBkgName_qqZZ = "{0}/Dbackground_qqZZ_{1}.root".format(opt.templateDir,appendName)
        bkgTempFile_qqZZ = ROOT.TFile(templateBkgName_qqZZ)
        bkgTemplate_qqZZ = bkgTempFile_qqZZ.Get("h_superDpsD")

        outTemplate_qqZZ = ROOT.TH1F("h_superDpsD","h_superDpsD",(bkgTemplate_qqZZ.GetNbinsX())*(bkgTemplate_qqZZ.GetNbinsY()),0,(bkgTemplate_qqZZ.GetNbinsX())*(bkgTemplate_qqZZ.GetNbinsY()))
        
        bin_qqZZ_tot=0
        for xbin in xrange(bkgTemplate_qqZZ.GetNbinsX()):
            for ybin in xrange(bkgTemplate_qqZZ.GetNbinsY()):
                outTemplate_qqZZ.SetBinContent(bin_qqZZ_tot+1,bkgTemplate_qqZZ.GetBinContent(xbin+1,ybin+1))
                bin_qqZZ_tot+=1
             

        OutTempBkgName_qqZZ = "{0}/Dbackground_qqZZ_{1}.root".format(opt.OutDir,appendName)
        OutTempFile_qqZZ = ROOT.TFile(OutTempBkgName_qqZZ,"RECREATE")
        outTemplate_qqZZ.Write()
        OutTempFile_qqZZ.Write()

        bkgTempFile_qqZZ.Close()
        OutTempFile_qqZZ.Close()

        templateBkgName_ggZZ = "{0}/Dbackground_ggZZ_{1}.root".format(opt.templateDir,appendName)
        bkgTempFile_ggZZ = ROOT.TFile(templateBkgName_ggZZ)
        bkgTemplate_ggZZ = bkgTempFile_ggZZ.Get("h_superDpsD")


        outTemplate_ggZZ = ROOT.TH1F("h_superDpsD","h_superDpsD",(bkgTemplate_ggZZ.GetNbinsX())*(bkgTemplate_ggZZ.GetNbinsY()),0,(bkgTemplate_ggZZ.GetNbinsX())*(bkgTemplate_ggZZ.GetNbinsY()))

        bin_ggZZ_tot=0
        for xbin in xrange(bkgTemplate_ggZZ.GetNbinsX()):
            for ybin in xrange(bkgTemplate_ggZZ.GetNbinsY()):
                outTemplate_ggZZ.SetBinContent(bin_ggZZ_tot+1,bkgTemplate_ggZZ.GetBinContent(xbin+1,ybin+1))
                bin_ggZZ_tot+=1

        OutTempBkgName_ggZZ = "{0}/Dbackground_ggZZ_{1}.root".format(opt.OutDir,appendName)
        OutTempFile_ggZZ = ROOT.TFile(OutTempBkgName_ggZZ,"RECREATE")
        outTemplate_ggZZ.Write()
        OutTempFile_ggZZ.Write()

        bkgTempFile_ggZZ.Close()
        OutTempFile_ggZZ.Close()

        templateBkgName_ZJetsCR_AllChans = "{0}/Dbackground_ZJetsCR_AllChans.root".format(opt.templateDir)
        bkgTempFile_ZJetsCR_AllChans = ROOT.TFile(templateBkgName_ZJetsCR_AllChans)
        bkgTemplate_ZJetsCR_AllChans = bkgTempFile_ZJetsCR_AllChans.Get("h_superDpsD")
       

        outTemplate_ZJetsCR_AllChans = ROOT.TH1F("h_superDpsD","h_superDpsD",(bkgTemplate_ZJetsCR_AllChans.GetNbinsX())*(bkgTemplate_ZJetsCR_AllChans.GetNbinsY()),0,(bkgTemplate_ZJetsCR_AllChans.GetNbinsX())*(bkgTemplate_ZJetsCR_AllChans.GetNbinsY()))
        
        bin_ZJetsCR_AllChans_tot=0
        for xbin in xrange(bkgTemplate_ZJetsCR_AllChans.GetNbinsX()):
            for ybin in xrange(bkgTemplate_ZJetsCR_AllChans.GetNbinsY()):
                outTemplate_ZJetsCR_AllChans.SetBinContent(bin_ZJetsCR_AllChans_tot+1,bkgTemplate_ZJetsCR_AllChans.GetBinContent(xbin+1,ybin+1))
                bin_ZJetsCR_AllChans_tot+=1

        OutTempBkgName_ZJetsCR_AllChans = "{0}/Dbackground_ZJetsCR_AllChans.root".format(opt.OutDir)
        OutTempFile_ZJetsCR_AllChans = ROOT.TFile(OutTempBkgName_ZJetsCR_AllChans,"RECREATE")
        outTemplate_ZJetsCR_AllChans.Write()
        OutTempFile_ZJetsCR_AllChans.Write()
        
        dataFileDir = "CMSdata"
        dataTreeName = "data_obs"
        dataFileName = "{0}/hzz{1}_{2}.root".format(dataFileDir,appendName,opt.lumi)
        data_obs_file = ROOT.TFile(dataFileName)
        data_obs_tree = data_obs_file.Get(dataTreeName)
        
        data_obs_2D = ROOT.TH2F("data_obs_2d","data_obs_2d",bkgTemplate_ZJetsCR_AllChans.GetNbinsX(),bkgTemplate_ZJetsCR_AllChans.GetXaxis().GetXbins().GetArray(),bkgTemplate_ZJetsCR_AllChans.GetNbinsY(),bkgTemplate_ZJetsCR_AllChans.GetYaxis().GetXbins().GetArray())
        
        for event in data_obs_tree:
            if (data_obs_tree.CMS_zz4l_mass < 141.0):
                if (data_obs_tree.CMS_zz4l_mass > 106.0):
                    data_obs_2D.Fill(data_obs_tree.CMS_zz4l_smd,data_obs_tree.CMS_zz4l_pseudoKD)

        data_obs_1D = ROOT.TH1F("h_superDpsD","h_superDpsD",(bin_tot),0,(bin_tot))
        bin_data_tot=0
        for xbin in xrange(data_obs_2D.GetNbinsX()):
            for ybin in xrange(data_obs_2D.GetNbinsY()):
                data_obs_1D.SetBinContent(bin_data_tot+1,data_obs_2D.GetBinContent(xbin+1,ybin+1))
                bin_data_tot+=1

        dataOutDir = "CMSdata1D"
        OutDataName = "{0}/hzz{1}_{2}.root".format(dataOutDir,appendName,opt.lumi)
        OutDataFile = ROOT.TFile(OutDataName,"RECREATE")
        data_obs_1D.Write()
        OutDataFile.Write()

        bkgTempFile_ZJetsCR_AllChans.Close()
        OutTempFile_ZJetsCR_AllChans.Close()
        data_obs_file.Close()
        OutDataFile.Close()
        
        c=c+1
    

    sys.exit()



# run the create_RM_cfg() as main()
if __name__ == "__main__":
    make1Dtemplates()


