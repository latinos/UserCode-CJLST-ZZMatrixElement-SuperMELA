#! /usr/bin/env python
import sys, os, pwd, commands
import optparse, shlex, re
import linecache
import math
from scipy.special import erf
from ROOT import *
import ROOT
from array import array
from systematicsClassSMD import *
from inputReader import *

##-------------------------------
## script to turn 1D templates and 2D cards/workspaces into 1D rootfiles and cards
##
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

    parser.add_option('-t', '--templateDir', type='string', dest='templateDir', default="templates2D" ,help="directory with 1D template histos")
    parser.add_option('-d', '--CardDir', type='string',dest='CardDir', default="cards_0-_8TeV/HCG/126/", help="directory with 2D cards and workspaces")
    parser.add_option('-o', '--OutDir', type='string', dest='OutDir', default='templates1D', help="directory to write 1D cards")
    parser.add_option('-l', '--Lumi', type='string', dest='lumi',default="19.63", help="lumi")
    parser.add_option('-e', '--Energy', type='string', dest='tev',default="8", help="sqrt(s) TeV")
    parser.add_option('-m', '--mH', type='float', dest='mH', default=126, help="CMS_zz4l_mass value")
    
    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()


# the main procedure
def make1DspinCards():
    
    # parse the arguments and options
    global opt, args
    parseOptions()

    ROOT.gSystem.AddIncludePath("-I$ROOFITSYS/include/")
    ROOT.gSystem.AddIncludePath("-Iinclude/")
    ROOT.gSystem.Load("libRooFit")
    ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
    
    c=1
    while (c < 4 ):
        appendName = ''
        if (c == 1): appendName = '4mu'
        elif (c == 2): appendName = '4e'
        elif (c == 3): appendName = '2e2mu'

        outName = "{0}/hzz4l_{1}S_{2}TeV.input.root".format(opt.OutDir,appendName,opt.tev)

        outFile = ROOT.TFile(outName,"RECREATE")
        
        inFileName = "{0}/hzz4l_{1}S_{2}TeV.input.root".format(opt.CardDir,appendName,opt.tev)

        inFile = ROOT.TFile(inFileName,"OPEN")

        w = inFile.Get("w")

        CMS_zz4l_mass = w.var("CMS_zz4l_mass")

        CMS_zz4l_mass.setVal(opt.mH)

        inDatacardName = "{0}/hzz4l_{1}S_{2}TeV_ALT.txt".format(opt.CardDir,appendName,opt.tev)
        outDatacardName = "{0}/hzz4l_{1}S_{2}TeV_ALT.txt".format(opt.OutDir,appendName,opt.tev)
        outDatacard = open(outDatacardName, "wb")

        normline = linecache.getline(inDatacardName,14)
        print normline

        ggH_norm = (float(normline.split()[1]))*(w.function("ggH_norm").getVal())
        print ggH_norm
        ggH_ALT_norm = (float(normline.split()[2]))*(w.function("ggH_ALT_norm").getVal())
        print ggH_ALT_norm
        bkg2d_qqzz_norm = float(normline.split()[3])
        print bkg2d_qqzz_norm
        bkg2d_ggzz_norm = float(normline.split()[4])
        print bkg2d_ggzz_norm
        bkg2d_zjets_norm = float(normline.split()[5])
        print bkg2d_zjets_norm

        for line in xrange(32):
            if line == 5:
                outDatacard.write("shapes * * hzz4l_{0}S_{1}TeV.input.root $PROCESS $PROCESS_$SYSTEMATIC \n".format(appendName,opt.tev))
                outDatacard.write("shapes data_obs a{0} hzz4l_{1}S_{2}TeV.input.root data_obs \n".format(c,appendName,opt.tev))
            elif line == 14:
                outDatacard.write("rate {0} {1} {2} {3} {4} \n".format(ggH_norm,ggH_ALT_norm,bkg2d_qqzz_norm,bkg2d_ggzz_norm,bkg2d_zjets_norm))
            elif line == 29 and c < 3:
                outDatacard.write("CMS_zz4l_smd_zjets_bkg_{0} shape - - - - 1.0 \n".format(c))
                outDatacard.write("CMS_zz4l_smd_leptScale_sig_{0} shape 1.0 1.0 - - - \n".format(c))
                outDatacard.write("CMS_zz4l_smd_leptResol_sig_{0} shape 1.0 1.0 - - - \n".format(c))
                break
            elif line == 30 and c == 3:
                outDatacard.write("CMS_zz4l_smd_zjets_bkg_{0} shape - - - - 1.0 \n".format(c))
                outDatacard.write("CMS_zz4l_smd_leptScale_sig_{0} shape 1.0 1.0 - - - \n".format(c))
                outDatacard.write("CMS_zz4l_smd_leptResol_sig_{0} shape 1.0 1.0 - - - \n".format(c))
                break
            else:
                outDatacard.write(linecache.getline(inDatacardName,line))
            line=line+1

        InSigName = "{0}/Dsignal_{1}.root".format(opt.templateDir,appendName)
        sigTempFile = ROOT.TFile(InSigName)
        sigTemplate = sigTempFile.Get("h_superDpsD")
        sigTemplate_syst1Up = sigTempFile.Get("h_superDpsD_LeptScaleUp")
        sigTemplate_syst1Down = sigTempFile.Get("h_superDpsD_LeptScaleDown")
        sigTemplate_syst2Up = sigTempFile.Get("h_superDpsD_LeptSmearUp")
        sigTemplate_syst2Down = sigTempFile.Get("h_superDpsD_LeptSmearDown")

        OutSigHist = sigTemplate.Clone('ggH')
        OutSigHist.SetNameTitle("ggH","ggH")
        OutSigHist.Scale(ggH_norm/OutSigHist.Integral("width"))
        outFile.cd()
        OutSigHist.Write()

        OutSigHist_syst1Up = sigTemplate_syst1Up.Clone('ggH_CMS_zz4l_smd_leptScale_sig_{0}Up'.format(c))
        OutSigHist_syst1Up.SetNameTitle('ggH_CMS_zz4l_smd_leptScale_sig_{0}Up'.format(c),'ggH_CMS_zz4l_smd_leptScale_sig_{0}Up'.format(c))
        OutSigHist_syst1Up.Scale(ggH_norm/OutSigHist_syst1Up.Integral("width"))
        outFile.cd()
        OutSigHist_syst1Up.Write()

        OutSigHist_syst1Down = sigTemplate_syst1Down.Clone('ggH_CMS_zz4l_smd_leptScale_sig_{0}Down'.format(c))
        OutSigHist_syst1Down.SetNameTitle('ggH_CMS_zz4l_smd_leptScale_sig_{0}Down'.format(c),'ggH_CMS_zz4l_smd_leptScale_sig_{0}Down'.format(c))
        OutSigHist_syst1Down.Scale(ggH_norm/OutSigHist_syst1Down.Integral("width"))
        outFile.cd()
        OutSigHist_syst1Down.Write()

        OutSigHist_syst2Up = sigTemplate_syst2Up.Clone('ggH_CMS_zz4l_smd_leptResol_sig_{0}Up'.format(c))
        OutSigHist_syst2Up.SetNameTitle('ggH_CMS_zz4l_smd_leptResol_sig_{0}Up'.format(c),'ggH_CMS_zz4l_smd_leptResol_sig_{0}Up'.format(c))
        OutSigHist_syst2Up.Scale(ggH_norm/OutSigHist_syst2Up.Integral("width"))
        outFile.cd()
        OutSigHist_syst2Up.Write()

        OutSigHist_syst2Down = sigTemplate_syst2Down.Clone('ggH_CMS_zz4l_smd_leptResol_sig_{0}Down'.format(c))
        OutSigHist_syst2Down.SetNameTitle('ggH_CMS_zz4l_smd_leptResol_sig_{0}Down'.format(c),'ggH_CMS_zz4l_smd_leptResol_sig_{0}Down'.format(c))
        OutSigHist_syst2Down.Scale(ggH_norm/OutSigHist_syst2Down.Integral("width"))
        outFile.cd()
        OutSigHist_syst2Down.Write()

        InSigName_ALT = "{0}/Dsignal_ALT_{1}.root".format(opt.templateDir,appendName)
        sigTempFile_ALT = ROOT.TFile(InSigName_ALT)
        sigTemplate_ALT = sigTempFile_ALT.Get("h_superDpsD")
        sigTemplate_ALT_syst1Up = sigTempFile_ALT.Get("h_superDpsD_LeptScaleUp")
        sigTemplate_ALT_syst1Down = sigTempFile_ALT.Get("h_superDpsD_LeptScaleDown")
        sigTemplate_ALT_syst2Up = sigTempFile_ALT.Get("h_superDpsD_LeptSmearUp")
        sigTemplate_ALT_syst2Down = sigTempFile_ALT.Get("h_superDpsD_LeptSmearDown")

        OutSigHist_ALT = sigTemplate_ALT.Clone('ggH_ALT')
        OutSigHist_ALT.SetNameTitle("ggH_ALT","ggH_ALT")
        OutSigHist_ALT.Scale(ggH_ALT_norm/OutSigHist_ALT.Integral("width"))
        outFile.cd()
        OutSigHist_ALT.Write()

        OutSigHist_ALT_syst1Up = sigTemplate_ALT_syst1Up.Clone('ggH_ALT_CMS_zz4l_smd_leptScale_sig_{0}Up'.format(c))
        OutSigHist_ALT_syst1Up.SetNameTitle('ggH_ALT_CMS_zz4l_smd_leptScale_sig_{0}Up'.format(c),'ggH_ALT_CMS_zz4l_smd_leptScale_sig_{0}Up'.format(c))
        OutSigHist_ALT_syst1Up.Scale(ggH_ALT_norm/OutSigHist_ALT_syst1Up.Integral("width"))
        outFile.cd()
        OutSigHist_ALT_syst1Up.Write()

        OutSigHist_ALT_syst1Down = sigTemplate_ALT_syst1Down.Clone('ggH_ALT_CMS_zz4l_smd_leptScale_sig_{0}Down'.format(c))
        OutSigHist_ALT_syst1Down.SetNameTitle('ggH_ALT_CMS_zz4l_smd_leptScale_sig_{0}Down'.format(c),'ggH_ALT_CMS_zz4l_smd_leptScale_sig_{0}Down'.format(c))
        OutSigHist_ALT_syst1Down.Scale(ggH_ALT_norm/OutSigHist_ALT_syst1Down.Integral("width"))
        outFile.cd()
        OutSigHist_ALT_syst1Down.Write()

        OutSigHist_ALT_syst2Up = sigTemplate_ALT_syst2Up.Clone('ggH_ALT_CMS_zz4l_smd_leptResol_sig_{0}Up'.format(c))
        OutSigHist_ALT_syst2Up.SetNameTitle('ggH_ALT_CMS_zz4l_smd_leptResol_sig_{0}Up'.format(c),'ggH_ALT_CMS_zz4l_smd_leptResol_sig_{0}Up'.format(c))
        OutSigHist_ALT_syst2Up.Scale(ggH_ALT_norm/OutSigHist_ALT_syst2Up.Integral("width"))
        outFile.cd()
        OutSigHist_ALT_syst2Up.Write()

        OutSigHist_ALT_syst2Down = sigTemplate_ALT_syst2Down.Clone('ggH_ALT_CMS_zz4l_smd_leptResol_sig_{0}Down'.format(c))
        OutSigHist_ALT_syst2Up.SetNameTitle('ggH_ALT_CMS_zz4l_smd_leptResol_sig_{0}Down'.format(c),'ggH_ALT_CMS_zz4l_smd_leptResol_sig_{0}Down'.format(c))
        OutSigHist_ALT_syst2Down.Scale(ggH_ALT_norm/OutSigHist_ALT_syst2Down.Integral("width"))
        outFile.cd()
        OutSigHist_ALT_syst2Down.Write()

        InBkgName_qqZZ = "{0}/Dbackground_qqZZ_{1}.root".format(opt.templateDir,appendName)
        bkgTempFile_qqZZ = ROOT.TFile(InBkgName_qqZZ)
        bkgTemplate_qqZZ = bkgTempFile_qqZZ.Get("h_superDpsD")

        OutBkgHist_qqZZ = bkgTemplate_qqZZ.Clone('bkg2d_qqzz')
        OutBkgHist_qqZZ.SetNameTitle("bkg2d_qqzz","bkg2d_qqzz")
        OutBkgHist_qqZZ.Scale(bkg2d_qqzz_norm/OutBkgHist_qqZZ.Integral("width"))
        outFile.cd()
        OutBkgHist_qqZZ.Write()

        InBkgName_ggZZ = "{0}/Dbackground_ggZZ_{1}.root".format(opt.templateDir,appendName)
        bkgTempFile_ggZZ = ROOT.TFile(InBkgName_ggZZ)
        bkgTemplate_ggZZ = bkgTempFile_ggZZ.Get("h_superDpsD")

        OutBkgHist_ggZZ = bkgTemplate_ggZZ.Clone('bkg2d_ggzz')
        OutBkgHist_ggZZ.SetNameTitle("bkg2d_ggzz","bkg2d_ggzz")
        OutBkgHist_ggZZ.Scale(bkg2d_ggzz_norm/OutBkgHist_ggZZ.Integral("width"))
        outFile.cd()
        OutBkgHist_ggZZ.Write()

        InBkgName_zjets = "{0}/Dbackground_ZJetsCR_AllChans.root".format(opt.templateDir)
        bkgTempFile_zjets = ROOT.TFile(InBkgName_zjets)
        
        bkgTemplate_zjets = bkgTempFile_qqZZ.Get("h_superDpsD")
        bkgTemplate_zjets_Up = bkgTempFile_qqZZ.Get("h_superDpsD")
        
        bkgTemplate_zjets_Down = bkgTempFile_zjets.Get("h_superDpsD")

        OutBkgHist_zjets = bkgTemplate_zjets.Clone('bkg2d_zjets')
        OutBkgHist_zjets.SetNameTitle("bkg2d_zjets","bkg2d_zjets")
        OutBkgHist_zjets.Scale(bkg2d_zjets_norm/OutBkgHist_zjets.Integral("width"))
        outFile.cd()
        OutBkgHist_zjets.Write()
        
        OutBkgHist_zjets_Up = bkgTemplate_zjets_Up.Clone('bkg2d_zjets_CMS_zz4l_smd_zjets_bkg_{0}Up'.format(c))
        OutBkgHist_zjets_Up.SetNameTitle('bkg2d_zjets_CMS_zz4l_smd_zjets_bkg_{0}Up'.format(c),'bkg2d_zjets_CMS_zz4l_smd_zjets_bkg_{0}Up'.format(c))
        OutBkgHist_zjets_Up.Scale(bkg2d_zjets_norm/OutBkgHist_zjets_Up.Integral("width"))
        outFile.cd()
        OutBkgHist_zjets_Up.Write()

        OutBkgHist_zjets_Down = bkgTemplate_zjets_Down.Clone('bkg2d_zjets_CMS_zz4l_smd_zjets_bkg_{0}Down'.format(c))
        OutBkgHist_zjets_Down.SetNameTitle('bkg2d_zjets_CMS_zz4l_smd_zjets_bkg_{0}Down'.format(c),'bkg2d_zjets_CMS_zz4l_smd_zjets_bkg_{0}Down'.format(c))
        OutBkgHist_zjets_Down.Scale(bkg2d_zjets_norm/OutBkgHist_zjets_Down.Integral("width"))
        outFile.cd()
        OutBkgHist_zjets_Down.Write()

        dataInDir = "CMSdata1D"
        InDataName = "{0}/hzz{1}_{2}.root".format(dataInDir,appendName,opt.lumi)
        InDataFile = ROOT.TFile(InDataName)
        InDataHist = InDataFile.Get("h_superDpsD")

        OutDataHist = InDataHist.Clone("data_obs")
        outFile.cd()
        OutDataHist.Write()

        inFile.Close()
        outFile.Close()
        outDatacard.close()

        c=c+1
    sys.exit()

# run the create_RM_cfg() as main()
if __name__ == "__main__":
    make1DspinCards()
