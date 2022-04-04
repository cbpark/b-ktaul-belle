#include <Math/Vector4D.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TTree.h>
#include <iostream>
#include <memory>  // std::shared_ptr
#include <tuple>   // std::tie
#include "input.h"
#include "variable.h"

using std::cout;

const double MINVISIBLE = 0.0;
const double PZTOT = 0.0;

const double EPSILON = 1.0e-10;
const unsigned int NEVAL = 5000;

int main(int argc, char *argv[]) {
    if (!(argc == 3)) {
        std::cerr
            << "usage: ./bin/mtau <event.root> <output.root>\n"
            << "  <event.root>: input root file (required).\n"
            << "  <output.root>: output file to store the result (required).\n";
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

    // the particle momenta from the input.
    Float_t px_ks, py_ks, pz_ks;
    Float_t px_mus, py_mus, pz_mus;
    Float_t px_htaus, py_htaus, pz_htaus;
    Float_t px_dt, py_dt, pz_dt;
    Float_t px_mut, py_mut, pz_mut;

    event->SetBranchAddress("px_ks", &px_ks);
    event->SetBranchAddress("py_ks", &py_ks);
    event->SetBranchAddress("pz_ks", &pz_ks);
    event->SetBranchAddress("px_mus", &px_mus);
    event->SetBranchAddress("py_mus", &py_mus);
    event->SetBranchAddress("pz_mus", &pz_mus);
    event->SetBranchAddress("px_htaus", &px_htaus);
    event->SetBranchAddress("py_htaus", &py_htaus);
    event->SetBranchAddress("pz_htaus", &pz_htaus);
    event->SetBranchAddress("px_dt", &px_dt);
    event->SetBranchAddress("py_dt", &py_dt);
    event->SetBranchAddress("pz_dt", &pz_dt);
    event->SetBranchAddress("px_mut", &px_mut);
    event->SetBranchAddress("py_mut", &py_mut);
    event->SetBranchAddress("pz_mut", &pz_mut);

    TFile outfile{argv[2], "recreate"};
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

    // the random number generator for mtau_random.
    auto rnd = std::make_shared<TRandom3>();
    rnd->SetSeed(42);

    const auto nentries = event->GetEntries();
    cout << "-- Total number of events: " << nentries << '\n';
#ifndef DEBUG
    // --------------------------------------------------------------------------
    // event loop
    for (auto iev = 0; iev != nentries; ++iev) {
#else
    auto itest = 7;
    for (auto iev = itest; iev != itest + 1; ++iev) {
#endif
        event->GetEntry(iev);
        // event->Show(iev);

        auto input = analysis::mkInputCM(
            {px_ks, py_ks, pz_ks}, {px_mus, py_mus, pz_mus},
            {px_htaus, py_htaus, pz_htaus}, {px_dt, py_dt, pz_dt},
            {px_mut, py_mut, pz_mut});

#ifdef DEBUG
        cout << "\nvis_sig: " << input.vis_sig() << '\n'
             << "vis_tag: " << input.vis_tag() << '\n'
             << "ptmiss: " << input.ptmiss() << "\n\n";
#endif

        cout << "---\n";
        cout << "let a1 = FourMomentum::new" << input.vis_sig() << ";\n";
        cout << "let a2 = FourMomentum::new" << input.vis_tag() << ";\n";
        cout << "let ptmiss = TransverseMomentum::new" << input.ptmiss()
             << ";\n";
        cout << "---\n";

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
        auto m2sb_sol = yam2::m2CConsSQP(input_kinematics, EPSILON, NEVAL);
        // input_kinematics.value().set_eps_constraint(1.0e-4);
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
    // event loop ends.
    // --------------------------------------------------------------------------

    infile.Close();

    // output->Print();
    output->Write();

    outfile.Close();
}
