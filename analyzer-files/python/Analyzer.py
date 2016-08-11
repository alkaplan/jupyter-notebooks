## This file is part of CMSOutreachExercise2011 derived from CMSOutreachExercise2010.
## Copyright (C) 2014 Instituto de Fisica de Cantabria and CERN.
## Based on the code of the CMSData Analysis School 2014 Long Exercise: 
## Search for the Higgs in ZZ -> 4 leptons decay channel (available 
## at https://github.com/bachtis/CMSDAS)

## CMSOutreachExercise2011 is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## CMSOutreachExercise2011 is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with CMSOutreachExercise2011. If not, see <http://www.gnu.org/licenses/>.

from DataFormats.FWLite import Events, Handle
from OutreachExercise2011.DecaysToLeptons import rootnotes
from ROOT import TFile, TTree
import ROOT
import os
import time
import array

def getFullPath(path):
    return os.path.join(os.environ['CMSSW_BASE'],
                        'src/OutreachExercise2011/DecaysToLeptons',
                        path)


class EventBox(object):
    def __init__(self):
        self.Z1_l1 = None
        self.Z1_l2 = None
        self.Z2_l1 = None
        self.Z2_l2 = None


class Object(object):
    """
    What is this?
    """
    def __init__(self, l1, l2=None):
        self.l1 = l1
        self.l2 = l2
        if l2:
            self.p4 = ROOT.TLorentzVector(self.l1.px() + self.l2.px(),
                                          self.l1.py() + self.l2.py(),
                                          self.l1.pz() + self.l2.pz(),
                                          self.l1.energy() + self.l2.energy())
        else:
            self.p4 = ROOT.TLorentzVector(self.l1.px(),
                                          self.l1.py(),
                                          self.l1.pz(),
                                          self.l1.energy())

    def mass(self):
        return self.p4.M()

    def pdgId(self):
        return 23

    def p4(self):
        return self.p4

    def pt(self):
        return self.p4.Pt()

    def px(self):
        return self.p4.Px()

    def py(self):
        return self.p4.Py()

    def pz(self):
        return self.p4.Pz()

    def energy(self):
        return self.p4.Energy()


class Analyzer (object):
    """
    Skeleton for building Analyzers
    """

    def __init__(self):
        self.vertexHandle = Handle('std::vector<reco::Vertex>')
        self.muonHandle = Handle('std::vector<pat::Muon>')
        self.electronHandle = Handle('std::vector<pat::Electron>')
        self.histograms = {}
        self.samples = ['data']
        self.plots = {}
        for sample in self.samples:
            self.histograms[sample] = {}
        self.data = []
        
        self.fileOut = TFile('output.root','recreate')
        self.t = TTree('t1','FourLepton')
        
        self.Lepton1_energy = array.array('f',[0])
        self.Lepton1_charge = array.array('f',[0])
        #self.Lepton1_global = array.array('i',[0])
        self.Lepton1_pt     = array.array('f',[0])
        self.Lepton1_px     = array.array('f',[0])
        self.Lepton1_py     = array.array('f',[0])
        self.Lepton1_pz     = array.array('f',[0])
        self.Lepton1_phi    = array.array('f',[0])
        self.Lepton1_eta    = array.array('f',[0])
        self.Lepton1_flavor = array.array('i',[0])
        
        self.Lepton2_energy = array.array('f',[0])
        self.Lepton2_charge = array.array('f',[0])
        #self.Lepton2_global = array.array('i',[0])
        self.Lepton2_pt     = array.array('f',[0])
        self.Lepton2_px     = array.array('f',[0])
        self.Lepton2_py     = array.array('f',[0])
        self.Lepton2_pz     = array.array('f',[0])
        self.Lepton2_phi    = array.array('f',[0])
        self.Lepton2_eta    = array.array('f',[0])
        self.Lepton2_flavor = array.array('i',[0])
        
        self.Lepton3_energy = array.array('f',[0])
        self.Lepton3_charge = array.array('f',[0])
        #self.Lepton3_global = array.array('i',[0])
        self.Lepton3_pt     = array.array('f',[0])
        self.Lepton3_px     = array.array('f',[0])
        self.Lepton3_py     = array.array('f',[0])
        self.Lepton3_pz     = array.array('f',[0])
        self.Lepton3_phi    = array.array('f',[0])
        self.Lepton3_eta    = array.array('f',[0])
        self.Lepton3_flavor = array.array('i',[0])
        
        self.Lepton4_energy = array.array('f',[0])
        self.Lepton4_charge = array.array('f',[0])
        #self.Lepton4_global = array.array('i',[0])
        self.Lepton4_pt     = array.array('f',[0])
        self.Lepton4_px     = array.array('f',[0])
        self.Lepton4_py     = array.array('f',[0])
        self.Lepton4_pz     = array.array('f',[0])
        self.Lepton4_phi    = array.array('f',[0])
        self.Lepton4_eta    = array.array('f',[0])
        self.Lepton4_flavor = array.array('i',[0])
        
        self.t.Branch('Lepton1_energy', self.Lepton1_energy, 'Lepton1_energy/F')
        self.t.Branch('Lepton1_charge', self.Lepton1_charge, 'Lepton1_charge/F')
        #self.t.Branch('Lepton1_global', self.Lepton1_global, 'Lepton1_global/F')
        self.t.Branch('Lepton1_pt',     self.Lepton1_pt,     'Lepton1_pt/F')
        self.t.Branch('Lepton1_px',     self.Lepton1_px,     'Lepton1_px/F')
        self.t.Branch('Lepton1_py',     self.Lepton1_py,     'Lepton1_py/F')
        self.t.Branch('Lepton1_pz',     self.Lepton1_pz,     'Lepton1_pz/F')
        self.t.Branch('Lepton1_phi',    self.Lepton1_phi,    'Lepton1_phi/F')
        self.t.Branch('Lepton1_eta',    self.Lepton1_eta,    'Lepton1_eta/F')
        self.t.Branch('Lepton1_flavor', self.Lepton1_flavor, 'Lepton1_flavor/I')
        
        self.t.Branch('Lepton2_energy', self.Lepton2_energy, 'Lepton2_energy/F')
        self.t.Branch('Lepton2_charge', self.Lepton2_charge, 'Lepton2_charge/F')
        #self.t.Branch('Lepton2_global', self.Lepton2_global, 'Lepton2_global/F')
        self.t.Branch('Lepton2_pt',     self.Lepton2_pt,     'Lepton2_pt/F')
        self.t.Branch('Lepton2_px',     self.Lepton2_px,     'Lepton2_px/F')
        self.t.Branch('Lepton2_py',     self.Lepton2_py,     'Lepton2_py/F')
        self.t.Branch('Lepton2_pz',     self.Lepton2_pz,     'Lepton2_pz/F')
        self.t.Branch('Lepton2_phi',    self.Lepton2_phi,    'Lepton2_phi/F')
        self.t.Branch('Lepton2_eta',    self.Lepton2_eta,    'Lepton2_eta/F')
        self.t.Branch('Lepton2_flavor', self.Lepton2_flavor, 'Lepton2_flavor/I')
        
        self.t.Branch('Lepton3_energy', self.Lepton3_energy, 'Lepton3_energy/F')
        self.t.Branch('Lepton3_charge', self.Lepton3_charge, 'Lepton3_charge/F')
        #self.t.Branch('Lepton3_global', self.Lepton3_global, 'Lepton3_global/F')
        self.t.Branch('Lepton3_pt',     self.Lepton3_pt,     'Lepton3_pt/F')
        self.t.Branch('Lepton3_px',     self.Lepton3_px,     'Lepton3_px/F')
        self.t.Branch('Lepton3_py',     self.Lepton3_py,     'Lepton3_py/F')
        self.t.Branch('Lepton3_pz',     self.Lepton3_pz,     'Lepton3_pz/F')
        self.t.Branch('Lepton3_phi',    self.Lepton3_phi,    'Lepton3_phi/F')
        self.t.Branch('Lepton3_eta',    self.Lepton3_eta,    'Lepton3_eta/F')
        self.t.Branch('Lepton3_flavor', self.Lepton3_flavor, 'Lepton3_flavor/I')
        
        self.t.Branch('Lepton4_energy', self.Lepton4_energy, 'Lepton4_energy/F')
        self.t.Branch('Lepton4_charge', self.Lepton4_charge, 'Lepton4_charge/F')
        #self.t.Branch('Lepton4_global', self.Lepton4_global, 'Lepton4_global/F')
        self.t.Branch('Lepton4_pt',     self.Lepton4_pt,     'Lepton4_pt/F')
        self.t.Branch('Lepton4_px',     self.Lepton4_px,     'Lepton4_px/F')
        self.t.Branch('Lepton4_py',     self.Lepton4_py,     'Lepton4_py/F')
        self.t.Branch('Lepton4_pz',     self.Lepton4_pz,     'Lepton4_pz/F')
        self.t.Branch('Lepton4_phi',    self.Lepton4_phi,    'Lepton4_phi/F')
        self.t.Branch('Lepton4_eta',    self.Lepton4_eta,    'Lepton4_eta/F')
        self.t.Branch('Lepton4_flavor', self.Lepton4_flavor, 'Lepton4_flavor/I')

    def muonID(self, muon, vertex):
        """
        Subclasses of Analyzer should define this method
        """
        return True

    def electronID(self, electron, vertex):
        """
        Subclasses of Analyzer should define this method
        """
        return True

    def leptonID(self, lepton, vertex):
        try:
            if abs(lepton.pdgId()) == 11:
                return self.electronID(lepton, vertex)
            elif abs(lepton.pdgId()) == 13:
                return self.muonID(lepton, vertex)
        except ZeroDivisionError:
            return False

    def readCollections(self, event, box, isFake=False):
        event.getByLabel('offlinePrimaryVertices', self.vertexHandle)
        event.getByLabel('patMuons', self.muonHandle)
        event.getByLabel('patElectrons', self.electronHandle)

        box.muons = self.muonHandle.product()
        box.electrons = self.electronHandle.product()
        box.vertex = self.vertexHandle.product()[0]

        #first select muons and electrons
        box.selectedMuons = []
        for mu in box.muons:
            if mu.pt() < 5 or abs(mu.eta()) > 2.4:
                continue
            try:
                if self.muonID(mu, box.vertex):
                    box.selectedMuons.append(mu)
            except ZeroDivisionError:
                continue

        box.selectedElectrons = []
        for ele in box.electrons:
            if ele.pt() < 5 or abs(ele.eta()) > 2.5:
                continue
            try:
                if self.electronID(ele, box.vertex):
                    box.selectedElectrons.append(ele)
            except ZeroDivisionError:
                continue

    def exportData(self):
        self.fileOut.Write()
        self.fileOut.Close()

    def analyze(self, box):
        return True

    def declareHistos(self):
        return True

    def fillHistos(self, box, sample, weight=1):
        return True

    def declareHisto(self, name, bins, min, max, xlabel=''):
        for sample in self.samples:
            self.histograms[sample][name] = ROOT.TH1F(sample + '_' + name,
                                                      name, bins, min, max)
            self.histograms[sample][name].GetXaxis().SetTitle(xlabel)
            self.histograms[sample][name].GetYaxis().SetTitle('events')

    def fillHisto(self, name, sample, value, weight=1):
        self.histograms[sample][name].Fill(value, weight)

    def addEvent(self, box):
        return

    def processSample(self, sample, maxEv=-1):
        print 'Processing Files'
        print sample.files

        events = Events(sample.files)
        print "%s events available for processing" % events.size()
        ts = time.time() 
        
        
        for N, event in enumerate(events):
            if N < 28451:
               continue
            if N > 28453 and N < 127159:
				continue
            if N > 127160 and N < 160330:
				continue
            if maxEv >= 0 and (N + 1) >= maxEv:
               break
            if N % 1000000 == 0:
                t2 = time.time()
                print "%s events processed in %s seconds" % (N + 1, t2 - ts)
            weight = 1
            box = EventBox()
            self.readCollections(event, box)
            if self.analyze(box) == False:
                continue
            #self.addEvent(box)
            t3 = time.time()
            print "Found one at %s after %s seconds!" % (N, t3 - ts)
            #self.fillHistos(box, sample.type, weight)
        tf = time.time()
        print "%s events processed in %s seconds" % (N + 1, tf - ts)

    def convertToPoisson(self, h):
        graph = ROOT.TGraphAsymmErrors()
        q = (1 - 0.6827) / 2.

        for i in range(1, h.GetNbinsX() + 1):
            x = h.GetXaxis().GetBinCenter(i)
            # XXX NOT USED!
            # xLow = h.GetXaxis().GetBinLowEdge(i)
            # xHigh = h.GetXaxis().GetBinUpEdge(i)
            y = h.GetBinContent(i)
            yLow = 0
            yHigh = 0
            if y != 0.0:
                yLow = y - ROOT.Math.chisquared_quantile_c(1 - q, 2 * y) / 2.
                yHigh = ROOT.Math.chisquared_quantile_c(q,
                                                        2 * (y + 1)) / 2. - y
                graph.SetPoint(i - 1, x, y)
                graph.SetPointEYlow(i - 1, yLow)
                graph.SetPointEYhigh(i - 1, yHigh)
                graph.SetPointEXlow(i - 1, 0.0)
                graph.SetPointEXhigh(i - 1, 0.0)

        graph.SetMarkerStyle(20)
        graph.SetLineWidth(2)
        graph.SetMarkerSize(1.)
        graph.SetMarkerColor(ROOT.kBlack)

        return graph

    def makePlot(self, histogram):
        # XXX NOT USED!
        # sandbox = []
        #canvas = ROOT.TCanvas(histogram)
        #if histogram in self.plots:
        #    return self.plots[histogram]['canvas']
        canvas = rootnotes.canvas(histogram)
        canvas.cd()
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptTitle(0)
        canvas.Range(-68.75, -7.5, 856.25, 42.5)
        canvas.SetFillColor(0)
        canvas.SetBorderMode(0)
        canvas.SetBorderSize(2)
        canvas.SetTickx(1)
        canvas.SetTicky(1)
        canvas.SetLeftMargin(0.15)
        canvas.SetRightMargin(0.05)
        canvas.SetTopMargin(0.05)
        canvas.SetBottomMargin(0.15)
        canvas.SetFrameFillStyle(0)
        canvas.SetFrameBorderMode(0)
        canvas.SetFrameFillStyle(0)
        canvas.SetFrameBorderMode(0)
        canvas.cd()

        frame = canvas.DrawFrame(
            self.histograms['data'][histogram].GetXaxis().GetXmin(),
            0.0, self.histograms['data'][histogram].GetXaxis().GetXmax(), 
            self.histograms['data'][histogram].GetMaximum()*1.2)

        frame.GetXaxis().SetLabelFont(42)
        frame.GetXaxis().SetLabelOffset(0.007)
        frame.GetXaxis().SetLabelSize(0.045)
        frame.GetXaxis().SetTitleSize(0.05)
        frame.GetXaxis().SetTitleOffset(1.15)
        frame.GetXaxis().SetTitleFont(42)
        frame.GetYaxis().SetLabelFont(42)
        frame.GetYaxis().SetLabelOffset(0.007)
        frame.GetYaxis().SetLabelSize(0.045)
        frame.GetYaxis().SetTitleSize(0.05)
        frame.GetYaxis().SetTitleOffset(1.4)
        frame.GetYaxis().SetTitleFont(42)
        frame.GetZaxis().SetLabelFont(42)
        frame.GetZaxis().SetLabelOffset(0.007)
        frame.GetZaxis().SetLabelSize(0.045)
        frame.GetZaxis().SetTitleSize(0.05)
        frame.GetZaxis().SetTitleFont(42)

        frame.GetXaxis().SetTitle(self.histograms['data'][histogram].GetXaxis().GetTitle())
        frame.GetYaxis().SetTitle(self.histograms['data'][histogram].GetYaxis().GetTitle())
        frame.Draw()

        self.histograms['data'][histogram].SetMarkerStyle(20)
        self.histograms['data'][histogram].Sumw2()
        dataG = None
        if self.histograms['data'][histogram].Integral() > 0:
            dataG = self.convertToPoisson(self.histograms['data'][histogram])
            dataG.Draw("Psame")

        legend = ROOT.TLegend(0.74, 0.84, 0.94, 0.94, "", "brNDC")
        legend.SetBorderSize(0)
        legend.SetLineColor(1)
        legend.SetLineStyle(1)
        legend.SetLineWidth(1)
        legend.SetFillColor(0)
        legend.SetFillStyle(0)
        legend.SetTextFont(42)
        legend.AddEntry(self.histograms['data'][histogram], "Data", "p")
        legend.Draw()

        pt = ROOT.TPaveText(0.1577181, 0.9562937, 0.9580537,
                            0.9947552, "brNDC")
        pt.SetBorderSize(0)
        pt.SetTextAlign(12)
        pt.SetFillStyle(0)
        pt.SetTextFont(42)
        pt.SetTextSize(0.03)
        pt.AddText(0.01, 0.4, "CMS")
        pt.Draw()

        self.plots[histogram] = {
            'canvas': canvas,
            'legend': legend,
            'dataG': dataG,
            'latex1': pt}

        canvas.RedrawAxis()
        canvas.Update()

        canvas.SaveAs(histogram+".png")          
       
        return canvas

    def makeAllPlots(self):
        return [self.makePlot(h) for h in self.histograms['data']]
