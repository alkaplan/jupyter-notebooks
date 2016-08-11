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

import itertools
import pprint
import array

from ROOT import *
from OutreachExercise2011.DecaysToLeptons.Analyzer import Analyzer, Object
from DataFormats.FWLite import Events, Handle

class FourLeptonAnalyzer(Analyzer):
    """
    XXX
    Some comment here.
    """
    def __init__(self):
        super(FourLeptonAnalyzer, self).__init__()

    #####CHANGE THIS METHOD TO CHANGE MUON ID######
    def muonID(self, muon, vertex):
        # if not (muon.isPFMuon() and (muon.isGlobalMuon() or
        #                              muon.isTrackerMuon())):
        # There are not PF muons on 2010 data! Only Global and Tracker muon
        if not (muon.isGlobalMuon() and muon.isTrackerMuon()):
            return False
        # higher-pt?
        if muon.pt() < 5 or abs(muon.eta()) > 2.4:
            return False
        # Position of the muon respect to the vertex
        if abs(muon.vertex().z() - vertex.z()) > 0.2:
            return False
        # OLD_Soft_Muon_selection_being_ph
        # The dB() method of the pat::Muon uses the version in IPTools,
        # so there are tiny differences between the values returned by
        # dxy(vertex->position()) and dB().
        # impact paramter dxy -> dB (w/ error edB) or ip (w/ error eip)
        # Impact parameter
        #if muon.innerTrack().dxy(vertex.position()) > 0.02:
        #    return False
        if muon.dB(muon.PV3D) > 0.02:
            return False

        # muon ISO variable
        #if (muon.chargedHadronIso() +
        #        max(0.0, muon.photonIso() + muon.neutralHadronIso() -
        #            0.5 * muon.puChargedHadronIso())) / muon.pt() > 0.6:
        # For 2010: I_trk + I_ECAL + I_HCAL
        # W&Z cross-section (2.9 pb-1) -> 0.15; CMSDAS 0.6
        if (muon.isolationR03().sumPt +
                muon.isolationR03().emEt +
                muon.isolationR03().hadEt) / muon.pt() > 0.15:
            return False

        # muon SIP variable # Symmetrized Impact Parameter in 2010?
        if (muon.dB(muon.PV3D) / muon.edB(muon.PV3D)) > 4:
            return False
        # chi2
        if muon.normChi2() > 10:
            return False
        # number of hits
        if muon.numberOfValidHits() < 10:
            return False

        return True

    #####CHANGE THIS METHOD TO CHANGE ELECTRON ID######
    def electronID(self, electron, vertex):
        return True
        if electron.pt() < 7 or abs(electron.eta()) > 2.5:
            return False

        #mvaRegions = [{'ptMin': 0, 'ptMax': 10, 'etaMin': 0.0,
        #               'etaMax': 0.8, 'mva': 0.47},
        #              {'ptMin': 0, 'ptMax': 10, 'etaMin': 0.8,
        #               'etaMax': 1.479, 'mva': 0.004},
        #              {'ptMin': 0, 'ptMax': 10, 'etaMin': 1.479,
        #               'etaMax': 3.0, 'mva': 0.295},
        #              {'ptMin': 10, 'ptMax': 99999999, 'etaMin': 0.0,
        #               'etaMax': 0.8, 'mva': -0.34},
        #              {'ptMin': 10, 'ptMax': 99999999, 'etaMin': 0.8,
        #               'etaMax': 1.479, 'mva': -0.65},
        #              {'ptMin': 10, 'ptMax': 99999999, 'etaMin': 1.479,
        #               'etaMax': 3.0, 'mva': 0.6},
        #              ]
        #ID = False
        #for region in mvaRegions:
        #    if electron.pt() >= region['ptMin'] and \
        #       electron.pt() < region['ptMax'] and \
        #       abs(electron.superCluster().eta()) >= region['etaMin'] and \
        #       abs(electron.superCluster().eta()) < region['etaMax'] and \
        #       electron.electronID("mvaNonTrigV0") > region['mva']:
        #       # "mvaNonTrigV0" not in 2010 data
        #        ID = True
        ##if not ID:

        #SimpleCutBasedEleId electron.electronID()
        # from top: 2010 - 2011 ???
        #0: fails
        #1: passes electron ID only
        #2: passes electron Isolation only
        #3: passes electron ID and Isolation only
        #4: passes conversion rejection
        #5: passes conversion rejection and ID
        #6: passes conversion rejection and Isolation
        #7: passes the whole selection

        if not electron.electronID("eidLoose"): 
            return False

        #photon conversion rejection? 

        if electron.gsfTrack().trackerExpectedHitsInner().numberOfHits() > 1:
            return False

        # electron ISO variable # For 2010: I_trk & I_ECAL & I_HCAL
        #if (electron.chargedHadronIso() + max(0.0,
        #        electron.photonIso() + electron.neutralHadronIso() -
        #        0.5 * electron.puChargedHadronIso())) / electron.pt() > 0.6:
        # barrel regions:
        if abs(electron.eta()) < 1.44:
            if electron.dr03TkSumPt() / electron.pt() > 0.09:
                return False
            if electron.dr03EcalRecHitSumEt() / electron.pt() > 0.07:
                return False
            if electron.dr03HcalTowerSumEt() / electron.pt() > 0.10:
                return False
        # endcap regions:
        if abs(electron.eta()) > 1.57:
            if electron.dr03TkSumPt() / electron.pt() > 0.04:
                return False
            if electron.dr03EcalRecHitSumEt() / electron.pt() > 0.05:
                return False
            if electron.dr03HcalTowerSumEt() / electron.pt() > 0.025:
                return False
        # excluded the region btw 1.44 and 1.57 (W&Z cross-section 2.9 pb-1)
        if abs(electron.eta()) > 1.44 and abs(electron.eta()) < 1.57:
                return False

        # similar to muons
        #if (electron.dr03TkSumPt() +
        #        electron.dr03EcalRecHitSumEt()  +
        #        electron.dr03HcalTowerSumEt()) / electron.pt() > 0.6:
        #    return False

        # electron SIP variable
        if (electron.dB(electron.PV3D) / electron.edB(electron.PV3D)) > 4:
            return False

        if electron.dB(electron.PV3D) > 0.04:
            return False

        # Position of the muon respect to the vertex
        if abs(electron.vertex().z() - vertex.z()) > 0.2:
            return False

        return True

    def select_zcandidates1(self, box):
        zcandidates = []
        for l1, l2 in itertools.combinations(box.leptons, 2):
            # they need to have same flavour and OS
            if abs(l1.pdgId()) != abs(l2.pdgId()):
                continue
            if l1.charge() + l2.charge() != 0:
                continue
            if (l1.pt() < 20 and l2.pt() < 10) or (l2.pt() < 20 and l1.pt() < 10):
                continue
            # now create a di lepton object and check mass
            z = Object(l1, l2)
            # Why this range: 12 -> 120? Ok first pair
            if not (z.mass() > 12 and z.mass() < 120):
                continue
            if z.mass() < 40:
                continue
            zcandidates.append(z)
        return zcandidates

    def select_zcandidates2(self, box):
        zcandidates = []
        for l1, l2 in itertools.combinations(box.leptons, 2):
            # they need to have same flavour and OS
            if abs(l1.pdgId()) != abs(l2.pdgId()):
                continue
            if l1.charge() + l2.charge() != 0:
                continue
            # now create a di lepton object and check mass
            z = Object(l1, l2)
            # Why this range: 12 -> 120? Not for second pair
            #if not (z.mass() > 12 and z.mass() < 120):
            #    continue
            if z.mass() < 4:
                continue
            zcandidates.append(z)
        return zcandidates

    #####ANALYSIS######
    def analyze(self, box):
        
        #print "I am working!"
        #####START FROM A bOX CONTAINING SELECTED MUONS AND ELECTRONS and MAKE
        #####FOUR LEPTON CANDIDATES

        # Now check if there are at least four leptons:
        box.leptons = set(box.selectedMuons + box.selectedElectrons)
        if len(box.leptons) < 4:
            return False

        # Now create Z candidates and apply cuts:
        box.zcandidates = self.select_zcandidates1(box)
        if len(box.zcandidates) == 0:
            return False

        # OK if there are more than one Z candidates
        # pick the one with the best mass
        sortedZs = sorted(box.zcandidates,
                          key=lambda x: abs(x.mass() - 91.118))
        box.Z1 = sortedZs[0]
        
        savedLeptons = [box.Z1.l1, box.Z1.l2]
        #print savedLeptons

        # now remove the used leptons from the list and make Z2 pairs
        box.leptons.remove(box.Z1.l1)
        box.leptons.remove(box.Z1.l2)

        # now the same thing with the second
        box.zcandidates2 = self.select_zcandidates2(box)
        if len(box.zcandidates2) == 0:
            return False

        # OK if there are more than one Z candidates
        # pick the one with the highest lepton pt sum
        sortedZ2s = sorted(box.zcandidates2,
                           key=lambda x: x.l1.pt() + x.l2.pt(), reverse=True)
        box.Z2 = sortedZ2s[0]
        
        savedLeptons.append(box.Z2.l1)
        savedLeptons.append(box.Z2.l2)
        
        self.addLeptonData(savedLeptons)

        # kill the candidate if a OS pair has mll<4 GeV
        for l1, l2 in itertools.combinations([box.Z1.l1, box.Z1.l2,
                                              box.Z2.l1, box.Z2.l2], 2):
            ll = Object(l1, l2)
            if (l1.charge() + l2.charge()) == 0:
                if ll.mass() < 4:
                    return False

        # create the ZZ
        box.ZZ = Object(box.Z1, box.Z2)
        return True

    def declareHistos(self):
        super(FourLeptonAnalyzer, self).declareHistos()

        ###ADD YOUR HISTOGRAMS AFTER THIS LINE AS AbOVE#####
        self.declareHisto('massZZ', 70, 0, 350, "m_{4l} [GeV]")
        self.declareHisto('massZ1', 20, 12, 120, "m_{Z1} [GeV]")
        self.declareHisto('massZ2', 20, 4, 74, "m_{Z2} [GeV]")

    def fillHistos(self, box, sample, weight=1):
        super(FourLeptonAnalyzer, self).fillHistos(box, sample, weight)

        self.fillHisto('massZZ', sample, box.ZZ.mass(), weight)
        self.fillHisto('massZ1', sample, box.ZZ.l1.mass(), weight)
        self.fillHisto('massZ2', sample, box.ZZ.l2.mass(), weight)

    def addEvent(self, box):
        self.data.append(box.ZZ.mass())
        #print 'working'
        #self.eventData.append(box)
        
    def addLeptonData(self, leptons):
    	l1 = leptons[0]
    	l2 = leptons[1]
    	l3 = leptons[2]
    	l4 = leptons[3]
    	    	
    	#lepton flavor declaration: if 1, it's a muon, if 2, it's an electron
    	for N, lept in enumerate(leptons):
			name = str(type(lept))[59:-2]
			if name == 'Electron':
				leptons[N].flavor = 2
			elif name == 'Muon':
				leptons[N].flavor = 1
				
    	
    	
    	self.Lepton1_energy[0] = l1.energy()
    	self.Lepton1_charge[0] = l1.charge()
       #self.Lepton1_global[0] = l1.isglobal()
    	self.Lepton1_pt[0]     = l1.pt()
    	self.Lepton1_px[0]     = l1.px()
    	self.Lepton1_py[0]     = l1.py()
    	self.Lepton1_pz[0]     = l1.pz()
    	self.Lepton1_phi[0]    = l1.phi()
    	self.Lepton1_eta[0]    = l1.eta()
    	self.Lepton1_flavor[0] = l1.flavor
    	
    	self.Lepton2_energy[0] = l2.energy()
    	self.Lepton2_charge[0] = l2.charge()
       #self.Lepton2_global[0] = l2.isglobal()
    	self.Lepton2_pt[0]     = l2.pt()
    	self.Lepton2_px[0]     = l2.px()
    	self.Lepton2_py[0]     = l2.py()
    	self.Lepton2_pz[0]     = l2.pz()
    	self.Lepton2_phi[0]    = l2.phi()
    	self.Lepton2_eta[0]    = l2.eta()
    	self.Lepton2_flavor[0] = l2.flavor
    	
    	self.Lepton3_energy[0] = l3.energy()
    	self.Lepton3_charge[0] = l3.charge()
       #self.Lepton3_global[0] = l3.isglobal()
    	self.Lepton3_pt[0]     = l3.pt()
    	self.Lepton3_px[0]     = l3.px()
    	self.Lepton3_py[0]     = l3.py()
    	self.Lepton3_pz[0]     = l3.pz()
    	self.Lepton3_phi[0]    = l3.phi()
    	self.Lepton3_eta[0]    = l3.eta()
    	self.Lepton3_flavor[0] = l3.flavor
    	
    	self.Lepton4_energy[0] = l4.energy()
    	self.Lepton4_charge[0] = l4.charge()
       #self.Lepton4_global[0] = l4.isglobal()
    	self.Lepton4_pt[0]     = l4.pt()
    	self.Lepton4_px[0]     = l4.px()
    	self.Lepton4_py[0]     = l4.py()
    	self.Lepton4_pz[0]     = l4.pz()
    	self.Lepton4_phi[0]    = l4.phi()
    	self.Lepton4_eta[0]    = l4.eta()
    	self.Lepton4_flavor[0] = l4.flavor
    	
    	self.t.Fill()
    	
    	


    	#self.eventData.append(event) '''
