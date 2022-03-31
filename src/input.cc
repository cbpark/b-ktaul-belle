#include "input.h"

#include <cmath>  // std::sqrt
#include <iostream>
#include "constant.h"

namespace analysis {
LorentzVector toLorentzVector(const Vector3 &v3, double mass) {
    double energy = std::sqrt(v3.mag2() + mass * mass);
    return {v3.x(), v3.y(), v3.z(), energy};
}

Input mkInputCM(const Vector3 &k_sig_v3, const Vector3 &mu_sig_v3,
                const Vector3 &htau_sig_v3, const Vector3 &d_tag_v3,
                const Vector3 &mu_tag_v3) {
    LorentzVector k_sig = toLorentzVector(k_sig_v3, MK);
    k_sig = boostToCM()(k_sig);

    LorentzVector mu_sig = toLorentzVector(mu_sig_v3, MMU);
    mu_sig = boostToCM()(mu_sig);

    LorentzVector htau_sig = toLorentzVector(htau_sig_v3, MPI);
    htau_sig = boostToCM()(htau_sig);

    LorentzVector d_tag = toLorentzVector(d_tag_v3, MD);
    d_tag = boostToCM()(d_tag);

    LorentzVector mu_tag = toLorentzVector(mu_tag_v3, MMU);
    mu_tag = boostToCM()(mu_tag);

    return {k_sig, mu_sig, htau_sig, d_tag, mu_tag};
}
}  // namespace analysis
