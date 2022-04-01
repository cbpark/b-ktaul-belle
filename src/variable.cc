#include "variable.h"

#include <TRandom3.h>
#include <YAM2/yam2.h>
#include <cmath>  // std::sqrt
#include <memory>
#include <optional>
#include "constant.h"
#include "input.h"

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

M2Reconstruction mkM2Reconstruction(
    const Input &input, const std::optional<yam2::M2Solution> &m2sol) {
    if (!m2sol) { return {-1.0, -1.0}; }

    auto m2sol_ = m2sol.value();
    auto k1sol = m2sol_.k1();
    LorentzVector k1sol_{k1sol.px(), k1sol.py(), k1sol.pz(), k1sol.e()};
    auto p_tau = input.htau_sig() + k1sol_;
    double mtau_sq = p_tau.mass2();

    // why negative?
    double mtau = mtau_sq < 0.0 ? -1.0 : std::sqrt(mtau_sq);

    return {m2sol_.m2(), mtau};
}

double mRecoilRandom(const Input &input, std::shared_ptr<TRandom3> rnd) {
    double cos_theta = rnd->Uniform(-1.0, 1.0);
    return mRecoil(input, cos_theta);
}
}  // namespace analysis
