# Compare hierarchy and validation QA histograms made
# by hierarchyPlots.py and validationPlots.py

import argparse
import array
import math
import os
import ROOT
import sys

class parameters(object):

    def __init__(self, hFileName, vFileName):
        # Hierarchy and validation histogram files
        self.hFileName = hFileName
        self.vFileName = vFileName


def compareHistos(pars):

    # Open MC histogram ROOT file
    print('plotMCHistos using {0} and {1}'.format(pars.hFileName, pars.vFileName))

    # Hierarchy histogram file
    hFile = ROOT.TFile.Open(pars.hFileName, 'read')
    # Validation histogram file
    vFile = ROOT.TFile.Open(pars.vFileName, 'read')

    # Define plotting canvas
    theCanvas = ROOT.TCanvas('theCanvas', '', 900, 700)
    ROOT.gROOT.SetStyle('Plain')
    ROOT.gStyle.SetOptStat(0)
    theCanvas.UseCurrentStyle()

    # For text labelling
    text = ROOT.TLatex()
    text.SetTextSize(0.075)
    text.SetNDC(True)

    # All interactions: muons, protons, piplus and piminus
    theCanvas.Divide(2,2)

    # Hits efficiency
    maxNHits = 1.0e4

    theCanvas.cd(1)
    h_muHitsEff = hFile.Get('muon_HitsEff')
    h_muHitsEff.GetXaxis().SetRangeUser(0.0, maxNHits)
    h_muHitsEff.GetYaxis().SetRangeUser(0.0, 1.01)
    h_muHitsEff.SetLineColor(ROOT.kBlue)
    h_muHitsEff.SetMarkerColor(ROOT.kBlue)
    ROOT.gPad.SetLogx()
    h_muHitsEff.Draw()
    v_muHitsEff = vFile.Get('muon_HitsEff')
    v_muHitsEff.SetLineColor(ROOT.kRed)
    v_muHitsEff.SetMarkerColor(ROOT.kRed)
    v_muHitsEff.Draw('same')
    text.DrawLatex(0.775, 0.25, '#mu')

    theCanvas.cd(2)
    h_pHitsEff = hFile.Get('proton_HitsEff')
    h_pHitsEff.GetXaxis().SetRangeUser(0.0, maxNHits)
    h_pHitsEff.GetYaxis().SetRangeUser(0.0, 1.01)
    h_pHitsEff.SetLineColor(ROOT.kBlue)
    h_pHitsEff.SetMarkerColor(ROOT.kBlue)
    ROOT.gPad.SetLogx()
    h_pHitsEff.Draw()
    v_pHitsEff = vFile.Get('proton_HitsEff')
    v_pHitsEff.SetLineColor(ROOT.kRed)
    v_pHitsEff.SetMarkerColor(ROOT.kRed)
    v_pHitsEff.Draw('same')
    text.DrawLatex(0.775, 0.25, 'p')

    theCanvas.cd(3)
    h_pipHitsEff = hFile.Get('piplus_HitsEff')
    h_pipHitsEff.GetXaxis().SetRangeUser(0.0, maxNHits)
    h_pipHitsEff.GetYaxis().SetRangeUser(0.0, 1.01)
    h_pipHitsEff.SetLineColor(ROOT.kBlue)
    h_pipHitsEff.SetMarkerColor(ROOT.kBlue)
    ROOT.gPad.SetLogx()
    h_pipHitsEff.Draw()
    v_pipHitsEff = vFile.Get('piplus_HitsEff')
    v_pipHitsEff.SetLineColor(ROOT.kRed)
    v_pipHitsEff.SetMarkerColor(ROOT.kRed)
    v_pipHitsEff.Draw('same')
    text.DrawLatex(0.775, 0.25, '#pi^{+}')

    theCanvas.cd(4)
    h_pimHitsEff = hFile.Get('piminus_HitsEff')
    h_pimHitsEff.GetXaxis().SetRangeUser(0.0, maxNHits)
    h_pimHitsEff.GetYaxis().SetRangeUser(0.0, 1.01)
    h_pimHitsEff.SetLineColor(ROOT.kBlue)
    h_pimHitsEff.SetMarkerColor(ROOT.kBlue)
    ROOT.gPad.SetLogx()
    h_pimHitsEff.Draw()
    v_pimHitsEff = vFile.Get('piminus_HitsEff')
    v_pimHitsEff.SetLineColor(ROOT.kRed)
    v_pimHitsEff.SetMarkerColor(ROOT.kRed)
    v_pimHitsEff.Draw('same')
    text.DrawLatex(0.775, 0.25, '#pi^{-}')

    theCanvas.Print('compare_allHitsEff.png')

    # Momentum efficiency
    maxMtm = 5.0

    theCanvas.cd(1)
    h_muMtmEff = hFile.Get('muon_MtmEff')
    h_muMtmEff.GetXaxis().SetRangeUser(0.0, maxMtm)
    h_muMtmEff.GetYaxis().SetRangeUser(0.0, 1.01)
    h_muMtmEff.SetLineColor(ROOT.kBlue)
    h_muMtmEff.SetMarkerColor(ROOT.kBlue)
    ROOT.gPad.SetLogx(0)
    h_muMtmEff.Draw()
    v_muMtmEff = vFile.Get('muon_MtmEff')
    v_muMtmEff.SetLineColor(ROOT.kRed)
    v_muMtmEff.SetMarkerColor(ROOT.kRed)
    v_muMtmEff.Draw('same')
    text.DrawLatex(0.775, 0.25, '#mu')

    theCanvas.cd(2)
    h_pMtmEff = hFile.Get('proton_MtmEff')
    h_pMtmEff.GetXaxis().SetRangeUser(0.0, maxMtm)
    h_pMtmEff.GetYaxis().SetRangeUser(0.0, 1.01)
    h_pMtmEff.SetLineColor(ROOT.kBlue)
    h_pMtmEff.SetMarkerColor(ROOT.kBlue)
    ROOT.gPad.SetLogx(0)
    h_pMtmEff.Draw()
    v_pMtmEff = vFile.Get('proton_MtmEff')
    v_pMtmEff.SetLineColor(ROOT.kRed)
    v_pMtmEff.SetMarkerColor(ROOT.kRed)
    v_pMtmEff.Draw('same')
    text.DrawLatex(0.775, 0.25, 'p')

    theCanvas.cd(3)
    h_pipMtmEff = hFile.Get('piplus_MtmEff')
    h_pipMtmEff.GetXaxis().SetRangeUser(0.0, maxMtm)
    h_pipMtmEff.GetYaxis().SetRangeUser(0.0, 1.01)
    h_pipMtmEff.SetLineColor(ROOT.kBlue)
    h_pipMtmEff.SetMarkerColor(ROOT.kBlue)
    ROOT.gPad.SetLogx(0)
    h_pipMtmEff.Draw()
    v_pipMtmEff = vFile.Get('piplus_MtmEff')
    v_pipMtmEff.SetLineColor(ROOT.kRed)
    v_pipMtmEff.SetMarkerColor(ROOT.kRed)
    v_pipMtmEff.Draw('same')
    text.DrawLatex(0.775, 0.25, '#pi^{+}')

    theCanvas.cd(4)
    h_pimMtmEff = hFile.Get('piminus_MtmEff')
    h_pimMtmEff.GetXaxis().SetRangeUser(0.0, maxMtm)
    h_pimMtmEff.GetYaxis().SetRangeUser(0.0, 1.01)
    h_pimMtmEff.SetLineColor(ROOT.kBlue)
    h_pimMtmEff.SetMarkerColor(ROOT.kBlue)
    ROOT.gPad.SetLogx(0)
    h_pimMtmEff.Draw()
    v_pimMtmEff = vFile.Get('piminus_MtmEff')
    v_pimMtmEff.SetLineColor(ROOT.kRed)
    v_pimMtmEff.SetMarkerColor(ROOT.kRed)
    v_pimMtmEff.Draw('same')
    text.DrawLatex(0.775, 0.25, '#pi^{-}')

    theCanvas.Print('compare_allMtmEff.png')

    # Close the histogram files
    hFile.Close()
    vFile.Close()


def run(args):

    pars = parameters(args.hFileName, args.vFileName)

    # Plot the comparisons
    compareHistos(pars)


def processArgs(parser):

    # Process script arguments
    parser.add_argument('--hFileName', default='MCHierarchy_Histos.root', metavar='fileName',
                        help='MC hierarchy histogram ROOT file [default "MCHierarchy_Histos.root"]')

    parser.add_argument('--vFileName', default='Validation_Histos.root', metavar='fileName',
                        help='Validation histogram ROOT file [default "Validation_Histos.root"]')


if __name__ == '__main__':

    # Process the command line arguments
    # Use "python hierarchyPlots.py --help" to see the full list
    parser = argparse.ArgumentParser(description='List of arguments')
    processArgs(parser)
    args = parser.parse_args()

    run(args)
