#include <Math/Vector3D.h>
#include <TFile.h>
#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TTree.h>
#include <iostream>
#include <memory>
#include "constant.h"
#include "input.h"
#include "variable.h"

using std::cout;
using Vector3 = ROOT::Math::XYZVector;

// a pair B mesons produced.
const double BB[2] = {analysis::MB, analysis::MB};
// signal side: B --> K mu tau
const double KMUTAU[3] = {analysis::MK, analysis::MMU, analysis::MTAU};
// tag side: B --> D mu nu
const double DMUNU[3] = {analysis::MD, analysis::MMU, 0.0};
// tau --> pi nu
const double PINU[2] = {analysis::MPI, 0.0};
// tau --> mu nu nu
const double MUNUNU[3] = {analysis::MMU, 0.0, 0.0};
// tau --> e nu nu
const double ENUNU[3] = {0.0, 0.0, 0.0};

const double MINVISIBLE = 0.0;
const double PZTOT = 0.0;

const double EPSILON = 1.0e-10;
const unsigned int NEVAL = 1000;

void print_momentum(const std::string pname, TLorentzVector &p) {
    cout << pname << ": ";
    p.Print();
}

int main(int, char *argv[]) {
    // total number of events
    const auto nev = 100000;

    TFile outfile{argv[1], "recreate"};
    cout << "-- The result will be stored in " << outfile.GetName() << '\n';

    // the event variables to calculate.
    Double_t mtau_random;
    Double_t m2s, mtau_m2s;
    Double_t m2sb, mtau_m2sb;

    auto output = std::make_shared<TTree>("var", "collider variables");
    output->Branch("mtau_random", &mtau_random, "mtau_random/D");
    output->Branch("m2s", &m2s, "m2s/D");
    output->Branch("mtau_m2s", &mtau_m2s, "mtau_m2s/D");
    output->Branch("m2sb", &m2sb, "m2sb/D");
    output->Branch("mtau_m2sb", &mtau_m2sb, "mtau_m2sb/D");

    // incoming electron and positron beams.
    TLorentzVector em{0.0, 0.0, analysis::EEM, analysis::EEM};
    TLorentzVector ep{0.0, 0.0, -analysis::EEP, analysis::EEP};
    TLorentzVector s = em + ep;

    // gRandom = new TRandom3(42);

    TGenPhaseSpace event;
    event.SetDecay(s, 2, BB);

    // the random number generator for mtau_random.
    auto rnd = std::make_shared<TRandom3>();
    rnd->SetSeed(42);

    for (auto iev = 0; iev != nev; ++iev) {
        // generate e+ e- --> B B event.
        event.Generate();
        auto b_sig = event.GetDecay(0);
        auto b_tag = event.GetDecay(1);

        TGenPhaseSpace b_sig_decay, b_tag_decay;
        // B(sig) --> K tau mu (three-body).
        b_sig_decay.SetDecay(*b_sig, 3, KMUTAU);
        // B(tag) --> D mu nu (three-body).
        b_tag_decay.SetDecay(*b_tag, 3, DMUNU);
        b_sig_decay.Generate();
        b_tag_decay.Generate();

        // the products of B(sig) decay.
        auto k_sig = b_sig_decay.GetDecay(0);
        auto mu_sig = b_sig_decay.GetDecay(1);
        auto tau_sig = b_sig_decay.GetDecay(2);

        // tau(sig) --> htau + nu.
        TGenPhaseSpace tau_sig_decay;
        // tau_sig_decay.SetDecay(*tau_sig, 2, PINU);
        // tau_sig_decay.SetDecay(*tau_sig, 3, MUNUNU);
        tau_sig_decay.SetDecay(*tau_sig, 3, ENUNU);
        tau_sig_decay.Generate();

        // the visible product of tau(sig) decay.
        auto htaus_sig = tau_sig_decay.GetDecay(0);

        // the products of B(tag) decay.
        auto d_tag = b_tag_decay.GetDecay(0);
        auto mu_tag = b_tag_decay.GetDecay(1);

        auto input = analysis::mkInputCM(
            {k_sig->Px(), k_sig->Py(), k_sig->Pz(), k_sig->E()},
            {mu_sig->Px(), mu_sig->Py(), mu_sig->Pz(), mu_sig->E()},
            {htaus_sig->Px(), htaus_sig->Py(), htaus_sig->Pz(), htaus_sig->E()},
            {d_tag->Px(), d_tag->Py(), d_tag->Pz(), d_tag->E()},
            {mu_tag->Px(), mu_tag->Py(), mu_tag->Pz(), mu_tag->E()});

        // mtau using random cos(theta).
        mtau_random = analysis::mRecoilRandom(input, rnd);

        // the input for calculating M2 variables (no vertex info).
        auto input_kinematics = input.to_input_kinematics(MINVISIBLE, PZTOT);

#ifdef DEBUG
        cout << "input kinematics:\n" << input_kinematics.value() << "\n\n";
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
        cout << m2sb_sol.value() << '\n';
#endif
        auto m2sb_rec = analysis::mkM2Reconstruction(input, m2sb_sol);
        std::tie(m2sb, mtau_m2sb) = m2sb_rec.get_result();

        // fill the event variables.
        output->Fill();
#ifdef DEBUG
        output->Show(0);
#endif

        // counter.
        auto nev = iev + 1;
        if ((nev < 10000 && nev % 1000 == 0) || nev % 10000 == 0) {
            cout << "... processed " << nev << " events.\n";
        }
    }

    // output->Print();
    output->Write();
}
