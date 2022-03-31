#include "input.h"

#include <YAM2/yam2.h>
#include <cmath>  // std::sqrt
#include <iostream>
#include <optional>
#include "constant.h"

namespace analysis {
yam2::FourMomentum fromLorentzVector(const LorentzVector &p) {
    return {p.e(), p.px(), p.py(), p.pz()};
}

yam2::TransverseMomentum fromVector2(const Vector2 &pt) {
    return {pt.x(), pt.y()};
}

std::optional<yam2::InputKinematics> Input::to_input_kinematics(
    double m_invisible, double pz_tot) const {
    auto p1 = fromLorentzVector(this->vis_sig());
    auto p2 = fromLorentzVector(this->vis_tag());

    // we have only one-step decay, so it's necessary to fake the second step.
    auto zero = yam2::FourMomentum();

    auto ptmiss = fromVector2(this->ptmiss());
    yam2::Mass m_parent{MB};

    return yam2::mkInput({p1, p2}, {zero, zero}, ptmiss,
                         yam2::Mass{m_invisible}, m_parent, {}, SQRTS,
                         {pz_tot});
}

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
