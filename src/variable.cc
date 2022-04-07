#include "variable.h"

#include <Math/Vector3D.h>
#include <TRandom3.h>
#include <YAM2/yam2.h>
#include <cmath>  // std::sqrt
#include <memory>
#include <optional>
#include "constant.h"
#include "input.h"

using Vector3 = ROOT::Math::XYZVector;

namespace analysis {
double mRecoil(const Input &input, double cos_theta) {
    auto kl_sig = input.k_sig() + input.mu_sig();
    double m_kl_sq = kl_sig.mass2();
    double e_kl = kl_sig.e();
    double p_kl = kl_sig.P();

    double m_tau_sq =
        MBSQ + m_kl_sq - 2.0 * (EBCM * e_kl + PBCM() * p_kl * cos_theta);
    return m_tau_sq < 0.0 ? -1.0 : std::sqrt(m_tau_sq);
}

double getCosTheta(const LorentzVector &p1, const LorentzVector &p2) {
    Vector3 p1_{p1.px(), p1.py(), p1.pz()};
    Vector3 p2_{p2.px(), p2.py(), p2.pz()};

    double norm1 = p1_.R();
    if (norm1 < 1.0e-10) { norm1 = 1.0e-10; }
    double norm2 = p2_.R();
    if (norm2 < 1.0e-10) { norm2 = 1.0e-10; }

    return p1_.Dot(p2_) / norm1 / norm2;
}

M2Reconstruction mkM2Reconstruction(
    const Input &input, const std::optional<yam2::M2Solution> &m2sol) {
    if (!m2sol) { return {-1.0, -1.0}; }

    auto m2sol_ = m2sol.value();

    /* reconstruct the tau mass: m_tau^2 = (htau + k1)^2. */
    // auto k1sol = m2sol_.k1();
    // LorentzVector k1sol_{k1sol.px(), k1sol.py(), k1sol.pz(), k1sol.e()};
    // auto p_tau = input.htau_sig() + k1sol_;
    // double mtau_sq = p_tau.mass2();
    // // why negative?
    // double mtau = mtau_sq < 0.0 ? -1.0 : std::sqrt(mtau_sq);

    /* reconstruct the tau mass using m_recoil. */
    auto k2sol = m2sol_.k2();
    // reconstructed invisible momentum on the tag side.
    LorentzVector k2sol_{k2sol.px(), k2sol.py(), k2sol.pz(), k2sol.e()};
    // reconstructed B_tag using the M2 solution.
    auto p_b_tag = input.vis_tag() + k2sol_;

    double cos_theta = getCosTheta(p_b_tag, input.kl_sig());
    double mtau = mRecoil(input, cos_theta);

    return {m2sol_.m2(), mtau};
}

double mRecoilRandom(const Input &input, std::shared_ptr<TRandom3> rnd) {
    double cos_theta = rnd->Uniform(-1.0, 1.0);
    return mRecoil(input, cos_theta);
}
}  // namespace analysis
