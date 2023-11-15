import awkward as ak
from pocket_coffea.lib.cut_definition import Cut

def dilepton(events, params, year, sample, **kwargs):
    MET = events[params["METbranch"][year]]
    # mass cut
    min_mass = events.ll.mass > 20
    #mass_window = events.ll.mass > 106 or events.ll.mass < 76
    mass_window = ak.any((events.ll.mass > 106) & (events.ll.mass < 76), axis=0)
    # Masks for same-flavor (SF) and opposite-sign (OS)
    SF = ((events.nMuonGood == 2) & (events.nElectronGood == 0)) | (
        (events.nMuonGood == 0) & (events.nElectronGood == 2)
    )
    OS = events.ll.charge == 0
    # SFOS = SF & OS
    not_SF = (events.nMuonGood == 1) & (events.nElectronGood == 1)
    mask = (
        (events.nLeptonGood == 2)
        & (ak.firsts(events.LeptonGood.pt) > params["pt_leading_lepton"])
        & (events.nJetGood >= params["njet"])
        & (events.nBJetGood >= params["nbjet"])
        & (MET.pt > params["met"])
        & OS
        & min_mass
        & mass_window
    )
    # Pad None values with False
    return ak.where(ak.is_none(mask), False, mask)

dilepton_presel = Cut(
    name="dilepton",
    params={
        "METbranch": {
            '2016_PreVFP': "MET",
            '2016_PostVFP': "MET",
            '2017': "MET",
            '2018': "MET",
        },
        "njet": 2,
        "nbjet": 0,
        "pt_leading_lepton": 15,
        "met": 40,
    },
    function=dilepton,
)
