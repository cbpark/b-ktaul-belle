#include <Math/Vector4D.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TTree.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>  // std::shared_ptr
#include <tuple>   // std::tie
#include "input.h"
#include "variable.h"

using std::cout;

const double MINVISIBLE = 0.0;
// const double PZTOT = 0.0;
const double PZTOT = analysis::PLONG;

const double SQRTS = 10.583;

const double EPSILON = 1.0e-8;
const unsigned int NEVAL = 1000;

int main(int argc, char *argv[]) {
    if (!(argc == 4)) {
        std::cerr
            << "usage: ./bin/mtau <event.root> <output.root> <tau_decay_mode>\n"
            << "  <event.root>: input root file (required).\n"
            << "  <output.root>: output file to store the result (required).\n"
            << "  <tau_decay_mode>: 0 (all) 1 (hadronic) 2 (leptonic)\n";
        return 1;
    }

    // the input root file.
    TFile infile{argv[1]};
    cout << "-- Input file: " << infile.GetName() << '\n';

    // check the trees in the input file.
    auto keys = infile.GetListOfKeys();
    if (keys->GetSize() < 1) {
        std::cerr << "-- The input file contains no tree.\n";
        infile.Close();
        return 1;
    }

    // get the tree.
    auto treename = infile.GetListOfKeys()->At(0)->GetName();
    cout << "-- The name of the input tree: " << treename << '\n';

    auto event = infile.Get<TTree>(treename);
    // event->Print();

    int tau_decay_mode = std::atoi(argv[3]);
    if (tau_decay_mode == 0) {
        cout << "-- Use all the decay modes of tau.\n";
    } else if (tau_decay_mode == 1) {
        cout << "-- Use hadronic decay modes of tau.\n";
    } else if (tau_decay_mode == 2) {
        cout << "-- Use leptonic decay modes of tau.\n";
    } else {
        std::cerr
            << "-- The tau decay mode must be set to be either 0, 1, or 2.\n";
        return 1;
    }

    // the particle momenta from the input.
    Float_t i_htaus;
    Float_t px_ks, py_ks, pz_ks, e_ks;
    Float_t px_mus, py_mus, pz_mus, e_mus;
    Float_t px_htaus, py_htaus, pz_htaus, e_htaus;
    Float_t px_dt, py_dt, pz_dt, e_dt;
    Float_t px_mut, py_mut, pz_mut, e_mut;

    // interaction point.
    Float_t vx_ip, vy_ip, vz_ip;

    // the secondary vertices
    Float_t vx_bs, vy_bs, vz_bs;
    Float_t vx_bt, vy_bt, vz_bt;

    cout << "-- Use the reco-level data.\n";
    event->SetBranchAddress("i_htaus", &i_htaus);
    event->SetBranchAddress("px_ks", &px_ks);
    event->SetBranchAddress("py_ks", &py_ks);
    event->SetBranchAddress("pz_ks", &pz_ks);
    event->SetBranchAddress("E_ks", &e_ks);
    event->SetBranchAddress("px_mus", &px_mus);
    event->SetBranchAddress("py_mus", &py_mus);
    event->SetBranchAddress("pz_mus", &pz_mus);
    event->SetBranchAddress("E_mus", &e_mus);
    event->SetBranchAddress("tpx_htau", &px_htaus);
    event->SetBranchAddress("tpy_htau", &py_htaus);
    event->SetBranchAddress("tpz_htau", &pz_htaus);
    event->SetBranchAddress("te_htaus", &e_htaus);
    event->SetBranchAddress("tpx_dt", &px_dt);
    event->SetBranchAddress("tpy_dt", &py_dt);
    event->SetBranchAddress("tpz_dt", &pz_dt);
    event->SetBranchAddress("te_dt", &e_dt);
    event->SetBranchAddress("px_mut", &px_mut);
    event->SetBranchAddress("py_mut", &py_mut);
    event->SetBranchAddress("pz_mut", &pz_mut);
    event->SetBranchAddress("E_mut", &e_mut);

    // event->SetBranchAddress("vx_ip", &vx_ip);
    // event->SetBranchAddress("vy_ip", &vy_ip);
    // event->SetBranchAddress("vz_ip", &vz_ip);
    // event->SetBranchAddress("vx_bs", &vx_bs);
    // event->SetBranchAddress("vy_bs", &vy_bs);
    // event->SetBranchAddress("vz_bs", &vz_bs);
    // event->SetBranchAddress("vx_bt", &vx_bt);
    // event->SetBranchAddress("vy_bt", &vy_bt);
    // event->SetBranchAddress("vz_bt", &vz_bt);

    event->SetBranchAddress("vx_ip", &vx_ip);
    event->SetBranchAddress("vy_ip", &vy_ip);
    event->SetBranchAddress("vz_ip", &vz_ip);
    event->SetBranchAddress("vx_bs", &vx_bs);
    event->SetBranchAddress("vy_bs", &vy_bs);
    event->SetBranchAddress("vz_bs", &vz_bs);
    event->SetBranchAddress("vx_bt", &vx_bt);
    event->SetBranchAddress("vy_bt", &vy_bt);
    event->SetBranchAddress("vz_bt", &vz_bt);

    TFile outfile{argv[2], "recreate"};
    cout << "-- The result will be stored in " << outfile.GetName() << '\n';

    // the event variables to calculate.
    Double_t mtau_random;
    Double_t m2s, mtau_m2s;
    Double_t m2sb, mtau_m2sb;
    Double_t m2sv_eq, mtau_m2sv_eq;
    Double_t m2v_eq, mtau_m2v_eq;

    auto output = std::make_shared<TTree>("var", "collider variables");
    output->Branch("mtau_random", &mtau_random, "mtau_random/D");
    output->Branch("m2s", &m2s, "m2s/D");
    output->Branch("mtau_m2s", &mtau_m2s, "mtau_m2s/D");
    output->Branch("m2sb", &m2sb, "m2sb/D");
    output->Branch("mtau_m2sb", &mtau_m2sb, "mtau_m2sb/D");
    output->Branch("m2sv_eq", &m2sv_eq, "m2sv_eq/D");
    output->Branch("mtau_m2sv_eq", &mtau_m2sv_eq, "mtau_m2sv_eq/D");
    output->Branch("m2v_eq", &m2v_eq, "m2v_eq/D");
    output->Branch("mtau_m2v_eq", &mtau_m2v_eq, "mtau_m2v_eq/D");

    LorentzVector k_sig, mu_sig, htau_sig;
    LorentzVector d_tag, mu_tag;
    Vector3 v_ip, v_bs, v_bt;

    // the random number generator for mtau_random.
    auto rnd = std::make_shared<TRandom3>();
    rnd->SetSeed(42);

    const auto nentries = event->GetEntries();
    cout << "-- Total number of events: " << nentries << '\n';
    unsigned int n_tau_decay_mode = 0;
#ifndef DEBUG
    // --------------------------------------------------------------------------
    // event loop
    for (auto iev = 0; iev != nentries; ++iev) {
#else
    auto itest = 0;
    for (auto iev = itest; iev != itest + 1; ++iev) {
#endif
        // counter.
        auto nev = iev + 1;
        if ((nev < 10000 && nev % 1000 == 0) || nev % 10000 == 0) {
            cout << "... processed " << nev << " events.\n";
        }

        event->GetEntry(iev);
        // event->Show(iev);

        if (tau_decay_mode != 0) {
            int id_tau_daughter = std::abs(i_htaus);
            if (id_tau_daughter == 211 || id_tau_daughter == 213 ||
                id_tau_daughter == 321) {
                if (tau_decay_mode != 1) { continue; }
            } else if (id_tau_daughter == 11 || id_tau_daughter == 13) {
                if (tau_decay_mode != 2) { continue; }
            } else {
                cout << "-- unknown i_htaus at event " << iev + 1 << ": "
                     << i_htaus << '\n';
                continue;
            }
        }
        ++n_tau_decay_mode;
#ifdef DEBUG
        cout << "-- id_tau_daughter: " << std::abs(i_htaus) << '\n';
#endif

        k_sig = LorentzVector(px_ks, py_ks, pz_ks, e_ks);
        mu_sig = LorentzVector(px_mus, py_mus, pz_mus, e_mus);
        htau_sig = LorentzVector(px_htaus, py_htaus, pz_htaus, e_htaus);
        d_tag = LorentzVector(px_dt, py_dt, pz_dt, e_dt);
        mu_tag = LorentzVector(px_mut, py_mut, pz_mut, e_mut);

        v_ip = Vector3(vx_ip, vy_ip, vz_ip);
        v_bs = Vector3(vx_bs, vy_bs, vz_bs);
        v_bt = Vector3(vx_bt, vy_bt, vz_bt);

        auto input = analysis::mkInput(k_sig, mu_sig, htau_sig, d_tag, mu_tag,
                                       {}, {v_ip}, {v_bs}, {v_bt});
        input.set_sqrt_s(SQRTS);
        input.set_pz_tot(PZTOT);

#ifdef DEBUG
        cout << "\nvis_sig: " << input.vis_sig() << '\n'
             << "vis_tag: " << input.vis_tag() << '\n'
             << "ptmiss: " << input.ptmiss() << "\n\n";

        // cout << "---\n";
        // cout << "let vis1 = FourMomentum::new" << input.vis_sig() << ";\n";
        // cout << "let vis2 = FourMomentum::new" << input.vis_tag() << ";\n";
        // cout << "let ptmiss = TransverseMomentum::new" << input.ptmiss()
        //      << ";\n";
        // cout << "let kl_sig = FourMomentum::new" << input.kl_sig() << ";\n";
        // cout << "---\n";
#endif

        // mtau using random cos(theta).
        mtau_random = analysis::mRecoilRandom(input, rnd);

        // the input for calculating M2 variables (no vertex info).
        auto input_kinematics = input.to_input_kinematics(MINVISIBLE);

#ifdef DEBUG
        // cout << "input kinematics:\n" << input_kinematics.value() << "\n\n";
#endif

        // reconstruction using M2s.
        auto m2s_sol = yam2::m2Cons(input_kinematics, EPSILON, NEVAL);
        auto m2s_rec = analysis::mkM2Reconstruction(input, m2s_sol);
        std::tie(m2s, mtau_m2s) = m2s_rec.get_result();

        // reconstruction using M2sB.
        auto m2sb_sol = yam2::m2CCons(input_kinematics, EPSILON, NEVAL);
        // input_kinematics.value().set_eps_constraint(1.0e-2);
        // auto m2sb_sol = yam2::m2CConsIneq(input_kinematics, EPSILON, NEVAL);
#ifdef DEBUG
        if (!m2sb_sol) {
            cout << "-- m2sb failed!\n";
        } else {
            cout << m2sb_sol.value() << '\n';
        }
#endif
        auto m2sb_rec = analysis::mkM2Reconstruction(input, m2sb_sol);
        std::tie(m2sb, mtau_m2sb) = m2sb_rec.get_result();

        auto input_kinematics_with_vertex =
            input.to_input_kinematics_with_vertex(input_kinematics, 0.0);
#ifdef DEBUG
        cout << "\ninput kinematics (with vertex):\n"
             << input_kinematics_with_vertex.value() << "\n\n";
#endif

        // reconstruction using M2sV(eq).
        // auto m2sv_eq_sol =
        //     yam2::m2ConsVertexEq(input_kinematics_with_vertex, EPSILON,
        //     NEVAL);
        auto m2sv_eq_sol =
            yam2::m2ConsVertexEq(input_kinematics_with_vertex, EPSILON, NEVAL);
        auto m2sv_eq_rec = analysis::mkM2Reconstruction(input, m2sv_eq_sol);
        std::tie(m2sv_eq, mtau_m2sv_eq) = m2sv_eq_rec.get_result();
#ifdef DEBUG
        if (!m2sv_eq_sol) {
            cout << "-- m2sv_eq failed!\n";
        } else {
            cout << m2sv_eq_sol.value() << '\n';
            cout << "mtot(m2sv_eq): " << analysis::mTotal(input, m2sv_eq_sol)
                 << "\n\n";
        }
#endif

        // reconstruction using M2sBV(eq).
        // auto m2sbv_eq_sol =
        //     yam2::m2CConsVertexEq(input_kinematics_with_vertex, EPSILON,
        //     NEVAL);
        // auto m2sbv_eq_rec = analysis::mkM2Reconstruction(input,
        // m2sbv_eq_sol); std::tie(m2sbv_eq, mtau_m2sbv_eq) =
        // m2sbv_eq_rec.get_result();

        // reconstruction using M2V(eq).
        auto m2v_eq_sol =
            yam2::m2VertexEq(input_kinematics_with_vertex, EPSILON, NEVAL);
        auto m2v_eq_rec = analysis::mkM2Reconstruction(input, m2v_eq_sol);
        std::tie(m2v_eq, mtau_m2v_eq) = m2v_eq_rec.get_result();
#ifdef DEBUG
        if (!m2v_eq_sol) {
            cout << "-- m2v_eq failed!\n";
        } else {
            cout << m2v_eq_sol.value() << '\n';
            cout << "mtot(m2v_eq): " << analysis::mTotal(input, m2v_eq_sol)
                 << '\n';
        }
#endif

        // fill the event variables.
        output->Fill();
#ifdef DEBUG
        // output->Show(0);
#endif
    }
    // event loop ends.
    // --------------------------------------------------------------------------

    infile.Close();

    // output->Print();
    output->Write();

    outfile.Close();

    cout << "-- " << n_tau_decay_mode << " events analyzed.\n";
}
