#include "input.h"

#include <cmath>
#include <iostream>
#include "constant.h"

namespace analysis {
LorentzVector toLorentzVector(const Vector3 &v3, double mass) {
    double energy = std::sqrt(v3.mag2() + mass * mass);
    return {v3.x(), v3.y(), v3.z(), energy};
}

Input mkInputCM(const Vector3 &k_sig_v3, const Vector3 &mu_sig_v3) {
    LorentzVector k_sig = toLorentzVector(k_sig_v3, MK);
    LorentzVector mu_sig = toLorentzVector(mu_sig_v3, MMU);

    return {k_sig, mu_sig};
}
}  // namespace analysis
