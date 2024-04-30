import awkward as ak
import sys
from pocket_coffea.workflows.base import BaseProcessorABC
from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.hist_manager import Axis
from pocket_coffea.lib.objects import (
    jet_correction,
    lepton_selection,
    jet_selection,
    btagging,
    bbtagging,
    get_dilepton,
    get_lepW,
    get_hadW
)
import numpy
import logging
from pocket_coffea.utils.logging import setup_logging
from pocket_coffea.lib.deltaR_matching import object_matching

class ttHbbBaseProcessor(BaseProcessorABC):
    def __init__(self, cfg: Configurator):
        super().__init__(cfg)

    def apply_object_preselection(self, variation):
        '''
        The ttHbb processor cleans
          - Electrons
          - Muons
          - Jets -> JetGood
          - BJet -> BJetGood

        '''
        # Include the supercluster pseudorapidity variable
        electron_etaSC = self.events.Electron.eta + self.events.Electron.deltaEtaSC
        self.events["Electron"] = ak.with_field(
            self.events.Electron, electron_etaSC, "etaSC"
        )
        # Build masks for selection of muons, electrons, jets, fatjets
        self.events["MuonGood"] = lepton_selection(self.events, "Muon", self.params)
        self.events["ElectronGood"] = lepton_selection(
            self.events, "Electron", self.params
        )
        leptons = ak.with_name(
            ak.concatenate((self.events.MuonGood, self.events.ElectronGood), axis=1),
            name='PtEtaPhiMCandidate',
        )
        self.events["LeptonGood"] = leptons[ak.argsort(leptons.pt, ascending=False)]


        self.events["FatJetGood"], self.jetGoodMask = jet_selection(
            self.events, "FatJet", self.params, "LeptonGood"
        )
        
        self.events["JetGood"], self.jetGoodMask = jet_selection(
            self.events, "Jet", self.params, "LeptonGood"
        )
        if self._isMC:
            self.events["GenJet"] = ak.with_field(self.events.GenJet, ak.zeros_like(self.events.GenJet.eta,dtype=int),"isFatJet")

            self.events["GenJetGood"], self.jetGoodMask = jet_selection(
            self.events, "GenJet", self.params
            )

            self.events["GenFatJet"] = ak.with_field(self.events.GenJetAK8, ak.ones_like(self.events.GenJetAK8.eta,dtype=int),"isFatJet")

            self.events["GenFatJetGood"], self.jetGoodMask = jet_selection(
            self.events, "GenFatJet", self.params
               )

        self.events["BJetGood"] = btagging(
            self.events["JetGood"], self.params.btagging.working_point[self._year]
        )

        isOld = False
        if ("v7" in self._dataset  ):
            isOld = True
        
        self.events["BBFatJetGoodT"] = bbtagging(
        self.events["FatJetGood"], self.params.btagging.working_point[self._year], "T", isOld
        )
        self.events["BBFatJetGoodM"] = bbtagging(
        self.events["FatJetGood"], self.params.btagging.working_point[self._year], "M", isOld
        )
        self.events["BBFatJetGoodL"] = bbtagging(
        self.events["FatJetGood"], self.params.btagging.working_point[self._year], "L", isOld
        )

        # to fix the errrors when running the dilepton MC
        #if len(self.events["FatJetGood"]) > 0 :
        #if self.params.btagging.working_point[self._year] in self.events["FatJetGood"]:
        #    if len(bbtagging( self.events["FatJetGood"], self.params.btagging.working_point[self._year], "T", isOld  )) > 0:
        #if "FatJet_"+self.params.btagging.working_point[self._year] in self.events :

        #else:
        #    self.events["BBFatJetGoodT"] =  None
        #    self.events["BBFatJetGoodM"] =  None
        #    self.events["BBFatJetGoodL"] =  None

        self.events["ll"] = get_dilepton(
            self.events.ElectronGood, self.events.MuonGood
        )
        self.events["lepW"] = get_lepW(
            self.events.ElectronGood, self.events.MuonGood, self.events.MET, self.events.JetGood, self.events.FatJetGood 
        )

        self.events["hadW"] = get_hadW(
            self.events.JetGood , self.events.BJetGood, self.events.FatJetGood
        )

    def count_objects(self, variation):
        self.events["nMuonGood"] = ak.num(self.events.MuonGood)
        self.events["nElectronGood"] = ak.num(self.events.ElectronGood)
        self.events["nLeptonGood"] = (
            self.events["nMuonGood"] + self.events["nElectronGood"]
        )
        self.events["nJetGood"] = ak.num(self.events.JetGood)
        self.events["nBJetGood"] = ak.num(self.events.BJetGood)
        self.events["nfatjet"]   = ak.num(self.events.FatJetGood)
        #if self.params.btagging.working_point[self._year] in self.events["FatJetGood"]:
        self.events["nBBFatJetGoodT"]   = ak.num(self.events.BBFatJetGoodT)
        self.events["nBBFatJetGoodM"]   = ak.num(self.events.BBFatJetGoodM)
        self.events["nBBFatJetGoodL"]   = ak.num(self.events.BBFatJetGoodL)
        #else:    
        #    self.events["nBBFatJetGoodT"]   = ak.num(0)
        #    self.events["nBBFatJetGoodM"]   = ak.num(0)
        #    self.events["nBBFatJetGoodL"]   = ak.num(0)       

    # Function that defines common variables employed in analyses and save them as attributes of `events`
    def define_common_variables_before_presel(self, variation):
        self.events["JetGood_Ht"] = ak.sum(abs(self.events.JetGood.pt), axis=1)
        if self._isMC:
            self.events["nJetGoodCFlavour"] = ak.sum(
                self.events.JetGood.hadronFlavour == 4, axis=1
            )
            self.events["nJetGoodBFlavour"] = ak.sum(
                self.events.JetGood.hadronFlavour == 5, axis=1
            ) 

    def process_extra_after_presel(self, variation):

        self.events["FatJetGood"] = ak.with_field(self.events.FatJetGood,2*numpy.log(self.events["FatJetGood"].msoftdrop/self.events["FatJetGood"].pt),"rhoQCD")
        self.events["BBFatJetGoodT"] = ak.with_field(self.events.BBFatJetGoodT,2*numpy.log(self.events["BBFatJetGoodT"].msoftdrop/self.events["BBFatJetGoodT"].pt),"rhoQCD")
        self.events["BBFatJetGoodM"] = ak.with_field(self.events.BBFatJetGoodM,2*numpy.log(self.events["BBFatJetGoodM"].msoftdrop/self.events["BBFatJetGoodM"].pt),"rhoQCD")
        self.events["BBFatJetGoodL"] = ak.with_field(self.events.BBFatJetGoodL,2*numpy.log(self.events["BBFatJetGoodL"].msoftdrop/self.events["BBFatJetGoodL"].pt),"rhoQCD")
        
        #logging.info("LHE: {}".format(self.events.LHEPart.pdgId[0]))
        #logging.info("LHE length: {}".format(len(self.events.LHEPart.pdgId)))
        '''
        isOutgoing = self.events.LHEPart.status == 1
        isParton = (abs(self.events.LHEPart.pdgId) < 6) | (
            self.events.LHEPart.pdgId == 21
        )
        quarks = self.events.LHEPart[isOutgoing & isParton]
        #logging.info("isOutgoing: {}".format(isOutgoing[0]))
        #logging.info("isParton: {}".format(isParton[0]))
        #logging.info("LHE outoing: {}".format(quarks.pdgId[0]))
        
        tops = self.events.GenPart[
            ((self.events.GenPart.pdgId == 6) | (self.events.GenPart.pdgId == -6))
            & (self.events.GenPart.hasFlags(['fromHardProcess']))
        ]
        #logging.info("tops {}".format(tops.pdgId[0]))
        onlytops = tops[tops.pdgId==6]
        genPartFromRad = onlytops.parent
        #logging.info("top parents {}".format(genPartFromRad.pdgId[0]))
        genPartFromRad = genPartFromRad.children
        genPartFromRad = genPartFromRad[((genPartFromRad.pdgId!=6) & (genPartFromRad.pdgId!=-6))]
        genPartFromRad = ak.flatten(genPartFromRad,axis=2)
        higgsProducts = genPartFromRad[genPartFromRad.pdgId==25]
        higgsProducts = higgsProducts.distinctChildrenDeep
        higgsProducts = ak.flatten(higgsProducts,axis=2)
        genPartFromRad = genPartFromRad[genPartFromRad.pdgId!=25]
        genPartFromRad = ak.concatenate((higgsProducts,genPartFromRad),axis=1)
        genPartFromRad = ak.with_field(genPartFromRad,0,"from_top")
        #logging.info("genPartFromRad {}".format(genPartFromRad.pdgId[0]))
        tops = tops[tops.hasFlags("isLastCopy")]
        children = tops.children
        #logging.info("children: {}".format(children.pdgId[0]))
        children = children[(children.pdgId!=6) & (children.pdgId!=-6)]
        children = ak.flatten(children,axis=2)
        genPartFromTop = children[(abs(children.pdgId)<6) | (children.pdgId==21)]
        #logging.info("genPartFromTop before loop {}".format(genPartFromTop.pdgId[0]))
        children = children[(abs(children.pdgId)>6) & (children.pdgId!=21) & ((abs(children.pdgId)<11) | (abs(children.pdgId)>18))]
        #logging.info("children before loop: {}".format(children.pdgId[0]))
        #while(len(ak.flatten(children,axis=1))):
        children = children.distinctChildrenDeep
        children = ak.flatten(children,axis=2)
        #logging.info("children in loop: {}".format(children.pdgId))
        genQuarks = children[(abs(children.pdgId)<6) | (children.pdgId==21)]
        genPartFromTop = ak.concatenate((genPartFromTop,genQuarks), axis=1)
        genPartFromTop = ak.with_field(genPartFromTop, 1, "from_top")
        #logging.info("genPartFromTop {}".format(genPartFromTop.pdgId))
        #children = children[(abs(children.pdgId)>6) & (children.pdgId!=21) & ((abs(children.pdgId)<11) | (abs(children.pdgId)>18))]

        genPart = ak.concatenate((genPartFromRad,genPartFromTop), axis=1)
        #logging.info("genjetgood {}".format(self.events.GenJetGood.eta))
        genJets = ak.concatenate((self.events.GenJetGood,self.events.GenFatJetGood), axis=1)
        #logging.info("genJets {}".format(genJets.eta))
        #logging.info("genJets isFatJet {}".format(genJets.isFatJet))
        #logging.info("genPart {}".format(genPart.pdgId[0]))
        #logging.info("genPart from top {}".format(genPart.from_top[0]))
        matched_jets, jetsBQuarks, jetsCQuarks, jetsFromTop = object_matching(
            genPart, genJets, dr_min=0.1
        )


        jetsBQuarks = ak.Array(jetsBQuarks)
        jetsCQuarks = ak.Array(jetsCQuarks)
        jetsFromTop = ak.Array(jetsFromTop)
        genJetsBQuarks = jetsBQuarks[genJets.isFatJet==0]
        genFatJetsBQuarks = jetsBQuarks[genJets.isFatJet==1]
        genJetsCQuarks = jetsCQuarks[genJets.isFatJet==0]
        genFatJetsCQuarks = jetsCQuarks[genJets.isFatJet==1]
        genJetsFromTop = jetsFromTop[genJets.isFatJet==0]
        genFatJetsFromTop = jetsFromTop[genJets.isFatJet==1]
        
        self.events["GenJetGood"] = ak.with_field(self.events.GenJetGood,genJetsBQuarks,"numBQuarksMatched")
        self.events["GenJetGood"] = ak.with_field(self.events.GenJetGood,genJetsCQuarks,"numCQuarksMatched")
        self.events["GenJetGood"] = ak.with_field(self.events.GenJetGood,genJetsFromTop,"isFromTop")
        self.events["GenFatJetGood"] = ak.with_field(self.events.GenFatJetGood,genFatJetsBQuarks,"numBQuarksMatched")
        self.events["GenFatJetGood"] = ak.with_field(self.events.GenFatJetGood,genFatJetsCQuarks,"numCQuarksMatched")
        self.events["GenFatJetGood"] = ak.with_field(self.events.GenFatJetGood,genFatJetsFromTop,"isFromTop")
        '''
        return

    def fill_histograms_extra(self, variation):
        '''
        This processor saves a metadata histogram with the number of
        events for chunk
        '''
        # # Filling the special histograms for events if they are present
        # if self._hasSubsamples:
        #     for subs in self._subsamples_names:
        #         if self._isMC & (
        #             "events_per_chunk" in self.hists_managers[subs].histograms
        #         ):
        #             hepc = self.hists_managers[subs].get_histogram("events_per_chunk")
        #             hepc.hist_obj.fill(
        #                 cat=hepc.only_categories[0],
        #                 variation="nominal",
        #                 year=self._year,
        #                 nEvents_initial=self.nEvents_initial,
        #                 nEvents_after_skim=self.nEvents_after_skim,
        #                 nEvents_after_presel=self.nEvents_after_presel,
        #             )
        #             self.output["processing_metadata"]["events_per_chunk"][
        #                 subs
        #             ] = hepc.hist_obj
        # else:
        #     if self._isMC & (
        #         "events_per_chunk" in self.hists_managers[self._sample].histograms
        #     ):
        #         hepc = self.hists_managers[self._sample].get_histogram(
        #             "events_per_chunk"
        #         )
        #         hepc.hist_obj.fill(
        #             cat=hepc.only_categories[0],
        #             variation="nominal",
        #             year=self._year,
        #             nEvents_initial=self.nEvents_initial,
        #             nEvents_after_skim=self.nEvents_after_skim,
        #             nEvents_after_presel=self.nEvents_after_presel,
        #         )
        #         self.output["processing_metadata"]["events_per_chunk"][
        #             self._sample
        #         ] = hepc.hist_obj
