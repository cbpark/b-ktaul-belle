#include <Math/Vector4D.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TTree.h>
#include <iostream>
#include <memory>
#include "input.h"
#include "variable.h"

using std::cout;

int main(int, char *argv[]) {
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

    Double_t m_tau_random;

    auto output = std::make_shared<TTree>("var", "collider variables");
    output->Branch("m_tau_random", &m_tau_random, "m_tau_random/D");

    auto rnd = std::make_shared<TRandom3>();
    rnd->SetSeed(42);

    const auto nentries = event->GetEntries();
    for (auto iev = 0; iev != nentries; ++iev) {
        event->GetEntry(iev);
        // event->Show(iev);

        auto input = analysis::mkInputCM(
            {px_ks, py_ks, pz_ks}, {px_mus, py_mus, pz_mus},
            {px_htaus, py_htaus, pz_htaus}, {px_dt, py_dt, pz_dt},
            {px_mut, py_mut, pz_mut});

        // cout << "k_sig: " << input.k_sig()
        //      << ", mass = " << input.k_sig().mass() << '\n';
        // cout << "mu_sig: " << input.mu_sig()
        //      << ", mass = " << input.mu_sig().mass() << '\n';
        // cout << "htau_sig: " << input.htau_sig()
        //      << ", mass = " << input.htau_sig().mass() << '\n';
        // cout << "d_tag: " << input.d_tag()
        //      << ", mass = " << input.d_tag().mass() << '\n';
        // cout << "mu_tag: " << input.mu_tag()
        //      << ", mass = " << input.mu_tag().mass() << '\n';
        // cout << "ptmiss: " << input.ptmiss() << '\n';

        m_tau_random = analysis::mRecoilRandom(input, rnd);
        // cout << "m_tau_random: " << m_tau_random << '\n';

        output->Fill();
    }

    infile.Close();
    cout << "-- Processed " << nentries << " events.\n";

    output->Print();
    output->Write();

    outfile.Close();
}
