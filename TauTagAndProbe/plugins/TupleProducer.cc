/*! Creates tuple for tau analysis.
This file is part of https://github.com/cms-tau-pog/TauTriggerTools. */

#include "Compression.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TauTriggerTools/Common/interface/CutTools.h"
#include "TauTriggerTools/Common/interface/GenTruthTools.h"
#include "TauTriggerTools/Common/interface/PatHelpers.h"
#include "TauTriggerTools/Common/interface/TriggerDescriptor.h"
#include "TauTriggerTools/TauTagAndProbe/interface/EventTuple.h"

namespace tau_trigger {

struct TupleProducerData {
    using Mutex = EventTuple::Mutex;
    using LockGuard = std::lock_guard<Mutex>;
    using SelectionHist = root_ext::SmartHistogram<cuts::ObjectSelector>;

    std::unique_ptr<EventTuple> eventTuple;
    std::unique_ptr<SelectionHist> selection;

    TupleProducerData(TFile& file)
    {
        eventTuple = std::make_unique<EventTuple>("events", &file, false);
        selection = std::make_unique<SelectionHist>("producer_selection");
        selection->SetOutputDirectory(&file);
    }
};

class TupleProducer : public edm::stream::EDProducer<edm::GlobalCache<TupleProducerData>> {
public:
    using SelectionHist = TupleProducerData::SelectionHist;
    using Cutter = cuts::Cutter<>;

    TupleProducer(const edm::ParameterSet& cfg, const TupleProducerData* producerData) :
        btagThreshold(cfg.getParameter<double>("btagThreshold")),
        isMC(cfg.getParameter<bool>("isMC")),
        triggerProcess(cfg.getParameter<std::string>("triggerProcess")),
        genEvent_token(mayConsume<GenEventInfoProduct>(cfg.getParameter<edm::InputTag>("genEvent"))),
        genParticles_token(mayConsume<std::vector<reco::GenParticle>>(cfg.getParameter<edm::InputTag>("genParticles"))),
        puInfo_token(mayConsume<std::vector<PileupSummaryInfo>>(cfg.getParameter<edm::InputTag>("puInfo"))),
        vertices_token(consumes<std::vector<reco::Vertex> >(cfg.getParameter<edm::InputTag>("vertices"))),
        signalMuon_token(consumes<pat::MuonRefVector>(cfg.getParameter<edm::InputTag>("signalMuon"))),
        taus_token(consumes<pat::TauCollection>(cfg.getParameter<edm::InputTag>("taus"))),
        jets_token(consumes<pat::JetCollection>(cfg.getParameter<edm::InputTag>("jets"))),
        met_token(consumes<pat::METCollection>(cfg.getParameter<edm::InputTag>("met"))),
        triggerResults_token(consumes<edm::TriggerResults>(edm::InputTag("TriggerResults", "", triggerProcess))),
        triggerDescriptors(cfg.getParameter<edm::VParameterSet>("hltPaths")),
        data(producerData),
        eventTuple(*data->eventTuple),
        selection(*data->selection)
    {
        produces<bool>();
    }

    static std::unique_ptr<TupleProducerData> initializeGlobalCache(const edm::ParameterSet&)
    {
        TFile& file = edm::Service<TFileService>()->file();
        file.SetCompressionAlgorithm(ROOT::kLZ4);
        file.SetCompressionLevel(4);
        return std::make_unique<TupleProducerData>(file);
    }

    static void globalEndJob(TupleProducerData* data)
    {
        TupleProducerData::LockGuard lock(data->eventTuple->GetMutex());
        data->eventTuple->Write();
        data->selection->WriteRootObject();
    }

private:
    static constexpr float default_value = ::tau_trigger::DefaultFillValue<float>();
    static constexpr int default_int_value = ::tau_trigger::DefaultFillValue<int>();
    static constexpr double deltaR2Thr = 0.5*0.5;

    virtual void beginRun(const edm::Run& run, const edm::EventSetup& setup)
    {
        HLTConfigProvider hltConfigProvider;
        bool changedConfig;
        if(!hltConfigProvider.init(run, setup, triggerProcess, changedConfig))
            throw exception("Unable to initialize HLTConfigProvider.");
        triggerDescriptors.updateGlobalIndices(hltConfigProvider.triggerNames());
    }

    virtual void produce(edm::Event& event, const edm::EventSetup&) override
    {
        event.put(std::make_unique<bool>(true));
        TupleProducerData::LockGuard lock(data->eventTuple->GetMutex());
        try {
            Cutter cut(&selection);
            fillTuple(event, cut);
        } catch(cuts::cut_failed&) {}
        selection.fill_selection();
    }

    void fillTuple(edm::Event& event, Cutter& cut)
    {
        cut(true, "total");
        eventTuple().run  = event.id().run();
        eventTuple().lumi = event.id().luminosityBlock();
        eventTuple().evt  = event.id().event();

        edm::Handle<std::vector<reco::Vertex>> vertices;
        event.getByToken(vertices_token, vertices);
        eventTuple().npv = static_cast<int>(vertices->size());

        edm::Handle<std::vector<reco::GenParticle>> hGenParticles;

        if(isMC) {
            edm::Handle<GenEventInfoProduct> genEvent;
            event.getByToken(genEvent_token, genEvent);
            eventTuple().genEventWeight = static_cast<float>(genEvent->weight());

            edm::Handle<std::vector<PileupSummaryInfo>> puInfo;
            event.getByToken(puInfo_token, puInfo);
            eventTuple().npu = gen_truth::GetNumberOfPileUpInteractions(puInfo);

            event.getByToken(genParticles_token, hGenParticles);
        }

        auto genParticles = hGenParticles.isValid() ? hGenParticles.product() : nullptr;
        std::vector<gen_truth::LeptonMatchResult> genLeptons;
        if(genParticles)
            genLeptons = gen_truth::CollectGenLeptons(*genParticles);

        edm::Handle<pat::MuonRefVector> signalMuonCollection;
        event.getByToken(signalMuon_token, signalMuonCollection);
        const pat::Muon* muon = signalMuonCollection.isValid() && !signalMuonCollection->empty()
                              ? &(*signalMuonCollection->at(0)) : nullptr;
        gen_truth::LeptonMatchResult gen_muon;
        LorentzVectorM muon_ref_p4;
        bool has_muon = false;
        if(muon) {
            gen_muon = gen_truth::LeptonGenMatch(muon->polarP4(), genLeptons);
            muon_ref_p4 = muon->polarP4();
            has_muon = true;
        } else {
            gen_muon = SelectGenLeg(genLeptons, false);
            if(gen_muon.match != GenLeptonMatch::NoMatch) {
                muon_ref_p4 = gen_muon.visible_p4;
                has_muon = true;
            }
        }
        cut(has_muon, "has_muon");

        // // Apply trigger match
        // edm::Handle<edm::TriggerResults> triggerResults;
        // event.getByToken(triggerResults_token, triggerResults);
        // const edm::TriggerNames& triggerNames = event.triggerNames(*triggerResults);

        edm::Handle<pat::METCollection> metCollection;
        event.getByToken(met_token, metCollection);
        const pat::MET& met = metCollection->at(0);
        const LorentzVectorM met_p4(met.pt(), 0, met.phi(), 0);
        eventTuple().met_pt = static_cast<float>(met.pt());
        eventTuple().met_phi = static_cast<float>(met.phi());
        eventTuple().muon_pt = muon ? static_cast<float>(muon->polarP4().pt()) : default_value;
        eventTuple().muon_eta = muon ? static_cast<float>(muon->polarP4().eta()) : default_value;
        eventTuple().muon_phi = muon ? static_cast<float>(muon->polarP4().phi()) : default_value;
        eventTuple().muon_mass = muon ? static_cast<float>(muon->polarP4().mass()) : default_value;
        eventTuple().muon_charge = muon ? muon->charge() : default_int_value;
        eventTuple().muon_iso = muon ? MuonIsolation(*muon) : default_value;
        eventTuple().muon_mt = muon ? Calculate_MT(muon->polarP4(), LorentzVectorM(met.pt(), 0, met.phi(), 0))
                                    : default_value;
        const bool has_gen_muon = gen_muon.match != GenLeptonMatch::NoMatch;
        eventTuple().muon_gen_match = static_cast<int>(gen_muon.match);
        eventTuple().muon_gen_charge = has_gen_muon ? gen_muon.gen_particle_lastCopy->charge() : default_int_value;
        eventTuple().muon_gen_vis_pt = has_gen_muon ? static_cast<float>(gen_muon.visible_p4.pt()) : default_value;
        eventTuple().muon_gen_vis_eta = has_gen_muon ? static_cast<float>(gen_muon.visible_p4.eta()) : default_value;
        eventTuple().muon_gen_vis_phi = has_gen_muon ? static_cast<float>(gen_muon.visible_p4.phi()) : default_value;
        eventTuple().muon_gen_vis_mass = has_gen_muon ? static_cast<float>(gen_muon.visible_p4.mass()) : default_value;

        edm::Handle<pat::TauCollection> taus;
        event.getByToken(taus_token, taus);

        edm::Handle<pat::JetCollection> jets;
        event.getByToken(jets_token, jets);

        const auto& selected_taus = CollectTaus(muon_ref_p4, *taus, genLeptons, deltaR2Thr);
        cut(!selected_taus.empty(), "has_tau");
        bool has_good_tau = false;
        for(const auto& tau_entry : selected_taus) {
            const pat::Tau* tau = tau_entry.reco_tau;
            const auto& gen_tau = tau_entry.gen_tau;
            const bool has_gen_tau = gen_tau.match != GenLeptonMatch::NoMatch;
            const LorentzVectorM tau_ref_p4 = tau ? tau->polarP4() : LorentzVectorM(gen_tau.visible_p4);
            if(!tau && !has_gen_tau)
                throw exception("Inconsistent tau entry");
            if(btagThreshold > 0  && !PassBtagVeto(muon_ref_p4, tau_ref_p4, *jets, btagThreshold, deltaR2Thr)) continue;

            eventTuple().tau_sel = tau_entry.selection;
            eventTuple().tau_pt = tau ? static_cast<float>(tau->polarP4().pt()) : default_value;
            eventTuple().tau_eta = tau ? static_cast<float>(tau->polarP4().eta()) : default_value;
            eventTuple().tau_phi = tau ? static_cast<float>(tau->polarP4().phi()) : default_value;
            eventTuple().tau_mass = tau ? static_cast<float>(tau->polarP4().mass()) : default_value;
            eventTuple().tau_charge = tau ? tau->charge() : default_int_value;

            eventTuple().tau_gen_match = static_cast<int>(gen_tau.match);
            eventTuple().tau_gen_charge = has_gen_tau ? gen_tau.gen_particle_firstCopy->charge() : default_int_value;
            eventTuple().tau_gen_vis_pt = has_gen_tau ? static_cast<float>(gen_tau.visible_p4.pt()) : default_value;
            eventTuple().tau_gen_vis_eta = has_gen_tau ? static_cast<float>(gen_tau.visible_p4.eta()) : default_value;
            eventTuple().tau_gen_vis_phi = has_gen_tau ? static_cast<float>(gen_tau.visible_p4.phi()) : default_value;
            eventTuple().tau_gen_vis_mass = has_gen_tau ? static_cast<float>(gen_tau.visible_p4.mass()) : default_value;
            eventTuple().tau_gen_rad_pt = has_gen_tau ? static_cast<float>(gen_tau.visible_rad_p4.pt()) : default_value;
            eventTuple().tau_gen_rad_eta = has_gen_tau ? static_cast<float>(gen_tau.visible_rad_p4.eta())
                                                       : default_value;
            eventTuple().tau_gen_rad_phi = has_gen_tau ? static_cast<float>(gen_tau.visible_rad_p4.phi())
                                                       : default_value;
            eventTuple().tau_gen_rad_energy = has_gen_tau ? static_cast<float>(gen_tau.visible_rad_p4.energy())
                                                        : default_value;
            eventTuple().tau_gen_n_charged_hadrons = has_gen_tau ? static_cast<int>(gen_tau.n_charged_hadrons)
                                                                 : default_int_value;
            eventTuple().tau_gen_n_neutral_hadrons = has_gen_tau ? static_cast<int>(gen_tau.n_neutral_hadrons)
                                                                 : default_int_value;
            eventTuple().tau_gen_n_gammas = has_gen_tau ? static_cast<int>(gen_tau.n_gammas) : default_int_value;
            eventTuple().tau_gen_n_gammas_rad = has_gen_tau ? static_cast<int>(gen_tau.n_gammas_rad)
                                                            : default_int_value;

            eventTuple().tau_decayMode = tau ? tau->decayMode() : default_int_value;
            eventTuple().tau_oldDecayModeFinding = tau ? tau->tauID("decayModeFinding") > 0.5f : default_int_value;

            for(const auto& tau_id_entry : tau_id::GetTauIdDescriptors()) {
                const auto& desc = tau_id_entry.second;
                desc.FillTuple(eventTuple, tau, default_value);
            }

            eventTuple().tau_dxy = tau ? tau->dxy() : default_value;
            eventTuple().tau_dxy_error = tau ? tau->dxy_error() : default_value;
            eventTuple().tau_ip3d = tau ? tau->ip3d() : default_value;
            eventTuple().tau_ip3d_error = tau ? tau->ip3d_error() : default_value;

            const pat::PackedCandidate* leadChargedHadrCand = tau
                    ? dynamic_cast<const pat::PackedCandidate*>(tau->leadChargedHadrCand().get())
                    : nullptr;
            eventTuple().tau_dz = leadChargedHadrCand ? leadChargedHadrCand->dz() : default_value;
            eventTuple().tau_dz_error = leadChargedHadrCand && leadChargedHadrCand->hasTrackDetails()
                    ? leadChargedHadrCand->dzError() : default_value;

            has_good_tau = true;
            eventTuple.Fill();
        }
        cut(has_good_tau, "btag_veto");
    }

private:
    const double btagThreshold;
    const bool isMC;
    const std::string triggerProcess;

    edm::EDGetTokenT<GenEventInfoProduct> genEvent_token;
    edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticles_token;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puInfo_token;
    edm::EDGetTokenT<std::vector<reco::Vertex>> vertices_token;
    edm::EDGetTokenT<pat::MuonRefVector> signalMuon_token;
    edm::EDGetTokenT<pat::TauCollection> taus_token;
    edm::EDGetTokenT<pat::JetCollection> jets_token;
    edm::EDGetTokenT<pat::METCollection> met_token;
    edm::EDGetTokenT<edm::TriggerResults> triggerResults_token;

    TriggerDescriptorCollection triggerDescriptors;
    const TupleProducerData* data;
    EventTuple& eventTuple;
    SelectionHist& selection;
};

} // namespace tau_trigger

#include "FWCore/Framework/interface/MakerMacros.h"
using TauTriggerTupleProducer = tau_trigger::TupleProducer;
DEFINE_FWK_MODULE(TauTriggerTupleProducer);
