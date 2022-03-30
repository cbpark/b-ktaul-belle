#include <Math/Vector4D.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include "input.h"

using std::cout;

int main(int, char *argv[]) {
    TFile infile{argv[1]};
#ifdef DEBUG
    cout << "Input file: " << infile.GetName() << '\n';
#endif  // DEBUG

    // check the trees in the input file.
    auto keys = infile.GetListOfKeys();
    if (keys->GetSize() < 1) {
        std::cerr << "The input file contains no tree.\n";
        infile.Close();
        return 1;
    }

    // get the tree.
    auto treename = infile.GetListOfKeys()->At(0)->GetName();
#ifdef DEBUG
    std::cout << "The name of the input tree: " << treename << '\n';
#endif

    auto event = infile.Get<TTree>(treename);
    // event->Print();

    Float_t px_ks, py_ks, pz_ks;
    Float_t px_mus, py_mus, pz_mus;
    Float_t px_htaus, py_htaus, pz_htaus;

    event->SetBranchAddress("px_ks", &px_ks);
    event->SetBranchAddress("py_ks", &py_ks);
    event->SetBranchAddress("pz_ks", &pz_ks);
    event->SetBranchAddress("px_mus", &px_mus);
    event->SetBranchAddress("py_mus", &py_mus);
    event->SetBranchAddress("pz_mus", &pz_mus);
    event->SetBranchAddress("px_htaus", &px_htaus);
    event->SetBranchAddress("py_htaus", &py_htaus);
    event->SetBranchAddress("pz_htaus", &pz_htaus);

    for (auto iev = 0; iev != 1; ++iev) {
        event->GetEntry(iev);
        // event->Show(iev);

        auto input = analysis::mkInputCM({px_ks, py_ks, pz_ks},
                                         {px_mus, py_mus, pz_mus});

        std::cout << "k_sig: " << input.k_sig()
                  << ", mass = " << input.k_sig().mass() << '\n';
        std::cout << "mu_sig: " << input.mu_sig()
                  << ", mass = " << input.mu_sig().mass() << '\n';
    }
}
