#include <Pythia8/Pythia.h>
#include <TFile.h>
#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TTree.h>
#include <cstdlib>
#include <iostream>
#include <vector>

using std::cout;

const double MK = 0.494;
const double MMU = 0.106;
const double MTAU = 1.777;
const double KMUTAU[3] = {MK, MMU, MTAU};

const double MPI = 0.140;
const double PINU[2] = {MPI, 0.0};
const double MUNUNU[3] = {MMU, 0.0, 0.0};

TLorentzVector getMomentum(const Pythia8::Particle &p) {
    return {p.px(), p.py(), p.pz(), p.e()};
}

int main(int argc, char *argv[]) {
    if (argc != 4) {
        std::cerr
            << "usage: ./gen_sig_bb <events.root> <nevent> <tau_decay>\n"
            << "  <events.root>: output file to store events.\n"
            << "  <nevent>: the number of events (not all events will be "
               "stored).\n"
            << "  <tau_decay>: 0 (tau -> pion nu) or 1 (tau -> l nu nu).\n";
        return 1;
    }

    Pythia8::Pythia pythia;

    pythia.readString("Random:setSeed = on");
#ifdef DEBUG
    pythia.readString("Random:seed = 41");
#else
    pythia.readString("Random:seed = 0");
#endif

    pythia.readString("Beams:frameType = 2");
    pythia.readString("Beams:idA =  11");
    pythia.readString("Beams:idB = -11");
    pythia.readString("Beams:eA = 8.0");
    pythia.readString("Beams:eB = 3.5");

    pythia.readString("PDF:lepton = off");
    // e+ e- --> gamma/Z (off-shell)
    pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
    pythia.readString("23:onMode = off");
    // gamma/Z --> b bbar
    pythia.readString("23:onIfAny = 5");
    // B+- --> D mu nu.
    pythia.readString("521:onMode = off");
    pythia.readString("521:onIfMatch = 14 -13 -421");

    pythia.init();

    // the root file to store events.
    TFile fout{argv[1], "recreate"};

    // the variables to store.
    int i_htaus;
    Float_t tpx_bs, tpy_bs, tpz_bs, te_bs;
    Float_t tpx_bt, tpy_bt, tpz_bt, te_bt;
    Float_t tpx_ks, tpy_ks, tpz_ks, te_ks;
    Float_t tpx_mus, tpy_mus, tpz_mus, te_mus;
    Float_t tpx_htau, tpy_htau, tpz_htau, te_htaus;
    Float_t tpx_dt, tpy_dt, tpz_dt, te_dt;
    Float_t tpx_mut, tpy_mut, tpz_mut, te_mut;

    auto tree = new TTree("h1", "");
    tree->Branch("i_htaus", &i_htaus, "i_htaus/I");
    tree->Branch("tpx_bs", &tpx_bs, "tpx_bs/F");
    tree->Branch("tpy_bs", &tpy_bs, "tpy_bs/F");
    tree->Branch("tpz_bs", &tpz_bs, "tpz_bs/F");
    tree->Branch("te_bs", &te_bs, "te_bs/F");
    tree->Branch("tpx_bt", &tpx_bt, "tpx_bt/F");
    tree->Branch("tpy_bt", &tpy_bt, "tpy_bt/F");
    tree->Branch("tpz_bt", &tpz_bt, "tpz_bt/F");
    tree->Branch("te_bt", &te_bt, "te_bt/F");
    tree->Branch("tpx_ks", &tpx_ks, "tpx_ks/F");
    tree->Branch("tpy_ks", &tpy_ks, "tpy_ks/F");
    tree->Branch("tpz_ks", &tpz_ks, "tpz_ks/F");
    tree->Branch("te_ks", &te_ks, "te_ks/F");
    tree->Branch("tpx_mus", &tpx_mus, "tpx_mus/F");
    tree->Branch("tpy_mus", &tpy_mus, "tpy_mus/F");
    tree->Branch("tpz_mus", &tpz_mus, "tpz_mus/F");
    tree->Branch("te_mus", &te_mus, "te_mus/F");
    tree->Branch("tpx_htau", &tpx_htau, "tpx_htau/F");
    tree->Branch("tpy_htau", &tpy_htau, "tpy_htau/F");
    tree->Branch("tpz_htau", &tpz_htau, "tpz_htau/F");
    tree->Branch("te_htaus", &te_htaus, "te_htaus/F");
    tree->Branch("tpx_dt", &tpx_dt, "tpx_dt/F");
    tree->Branch("tpy_dt", &tpy_dt, "tpy_dt/F");
    tree->Branch("tpz_dt", &tpz_dt, "tpz_dt/F");
    tree->Branch("te_dt", &te_dt, "te_dt/F");
    tree->Branch("tpx_mut", &tpx_mut, "tpx_mut/F");
    tree->Branch("tpy_mut", &tpy_mut, "tpy_mut/F");
    tree->Branch("tpz_mut", &tpz_mut, "tpz_mut/F");
    tree->Branch("te_mut", &te_mut, "te_mut/F");

    TLorentzVector b_sig, b_tag;
    TLorentzVector *k_sig, *mu_sig, *tau_sig, *htau_sig;
    TLorentzVector d_tag, mu_tag;

#ifdef DEBUG
    gRandom = new TRandom3(42);
#endif

    // Total number of events to generate.
    //
    // note that the actual number of events could be different due to
    // e+ e- --> B0 B0 events, which we will skip.
    int nevent = std::atoi(argv[2]);
    int igen = 0;
    int tau_mode = std::atoi(argv[3]);
    if (tau_mode != 0 && tau_mode != 1) {
        std::cerr << "tau_decay must be either 0 or 1.\n";
        return 1;
    }
    // --------------------------------------------------------------------------
    // event loop
    for (int iev = 0; iev < nevent; ++iev) {
        if (!pythia.next()) { continue; }

        auto event = pythia.event;
        std::vector<int> bmeson_list;
        for (auto iptl = 0; iptl != event.size(); ++iptl) {
            auto particle = event[iptl];
            // search for B+ B- mesons.
            if (particle.statusAbs() == 82 && particle.idAbs() == 521) {
                bmeson_list.push_back(iptl);
            }
        }

        // skip if the event doesn't have a B-meson pair
        if (bmeson_list.size() != 2) { continue; }

        // skip unless B+ --> D mu nu.
        auto b_tag_daughters_list = event[bmeson_list[0]].daughterList();
        if (b_tag_daughters_list.size() != 3) { continue; }

        b_tag = getMomentum(event[bmeson_list[0]]);
        auto idgt_ok = 0;
        for (const auto &idgt : b_tag_daughters_list) {
            auto b_tag_daughter = event[idgt];
            if (b_tag_daughter.idAbs() == 13) {
                // muon(tag).
                mu_tag = getMomentum(b_tag_daughter);
                ++idgt_ok;
            } else if (b_tag_daughter.idAbs() == 421) {
                // D(tag).
                d_tag = getMomentum(b_tag_daughter);
                ++idgt_ok;
            }
        }

        // skip if the decay products of B(tag) don't have a muon and a D meson.
        if (idgt_ok != 2) { continue; }

        // B(sig).
        b_sig = getMomentum(event[bmeson_list[1]]);

        // decay B(sig) using TGenPhaseSpace.
        TGenPhaseSpace b_sig_decay;
        // B(sig) --> K mu tau (three-body).
        b_sig_decay.SetDecay(b_sig, 3, KMUTAU);
        b_sig_decay.Generate();

        k_sig = b_sig_decay.GetDecay(0);
        mu_sig = b_sig_decay.GetDecay(1);
        tau_sig = b_sig_decay.GetDecay(2);

        TGenPhaseSpace tau_sig_decay;
        if (tau_mode == 0) {
            // tau(sig) --> pi nu
            tau_sig_decay.SetDecay(*tau_sig, 2, PINU);
        } else {
            // tau(sig) --> mu nu nu
            tau_sig_decay.SetDecay(*tau_sig, 3, MUNUNU);
        }
        tau_sig_decay.Generate();

        htau_sig = tau_sig_decay.GetDecay(0);

        if (tau_mode == 0) {
            i_htaus = 211;  // pion
        } else {
            i_htaus = 13;  // muon
        }
        tpx_bs = b_sig.Px();
        tpy_bs = b_sig.Py();
        tpz_bs = b_sig.Pz();
        te_bs = b_sig.E();
        tpx_bt = b_tag.Px();
        tpy_bt = b_tag.Py();
        tpz_bt = b_tag.Pz();
        te_bt = b_tag.E();
        tpx_ks = k_sig->Px();
        tpy_ks = k_sig->Py();
        tpz_ks = k_sig->Pz();
        te_ks = k_sig->E();
        tpx_mus = mu_sig->Px();
        tpy_mus = mu_sig->Py();
        tpz_mus = mu_sig->Pz();
        te_mus = mu_sig->E();
        tpx_htau = htau_sig->Px();
        tpy_htau = htau_sig->Py();
        tpz_htau = htau_sig->Pz();
        te_htaus = htau_sig->E();
        tpx_dt = d_tag.Px();
        tpy_dt = d_tag.Py();
        tpz_dt = d_tag.Pz();
        te_dt = d_tag.E();
        tpx_mut = mu_tag.Px();
        tpy_mut = mu_tag.Py();
        tpz_mut = mu_tag.Pz();
        te_mut = mu_tag.E();

        // fill in the entries.
        tree->Fill();

        ++igen;
    }
    // event loop ends.
    // --------------------------------------------------------------------------

    tree->Write();
    fout.Close();

    pythia.stat();

    cout << "-- done.\n";
    cout << "-- the number of events: " << igen << '\n';
    cout << "-- the generated events have been stored in " << fout.GetName()
         << '\n';
}
