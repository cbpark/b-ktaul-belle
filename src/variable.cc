#include "variable.h"

#include <TRandom3.h>
#include <YAM2/yam2.h>
#include <cmath>  // std::sqrt
#include <memory>
#include <optional>
#include "constant.h"
#include "input.h"

namespace analysis {
M2Reconstruction mkM2Reconstruction(
    const Input &input, const std::optional<yam2::M2Solution> &m2sol) {
    if (!m2sol) { return {-1.0, -1.0}; }

    auto m2sol_ = m2sol.value();
    auto k2sol = m2sol_.k2();
    // reconstructed invisible momentum on the tag side.
    LorentzVector k2sol_{k2sol.px(), k2sol.py(), k2sol.pz(), k2sol.e()};

    // reconstructed B_tag using the M2 solution.
    auto p_b_tag = input.vis_tag() + k2sol_;

    // p(CM) = p(B_sig) + p(B_tag)
    //       = p(K_sig l_sig) + p(tau_sig) + p(B_tag).
    // therefore, p(tau_sig) = p(CM) - p(K_sig l_sig) - p(B_tag).
    auto p_tau = CM() - input.kl_sig() - p_b_tag;
    double mtau_sq = p_tau.mass2();

    // why negative?
    double mtau = mtau_sq < 0.0 ? -1.0 : std::sqrt(mtau_sq);

    return {m2sol_.m2(), mtau};
}

double mRecoilRandom(const Input &input, std::shared_ptr<TRandom3> rnd) {
    auto kl_sig = input.k_sig() + input.mu_sig();
    double m_kl_sq = kl_sig.mass2();
    double e_kl = kl_sig.e();
    double p_kl = kl_sig.P();

    double m_b_sq = MB * MB;
    double e_b = 0.5 * SQRTS;
    double p_b = std::sqrt(e_b * e_b - m_b_sq);

    double cos_theta = rnd->Uniform(-1.0, 1.0);

    double m_tau_sq =
        m_b_sq + m_kl_sq - 2.0 * (e_b * e_kl + p_b * p_kl * cos_theta);
    return m_tau_sq < 0.0 ? -1.0 : std::sqrt(m_tau_sq);
}
}  // namespace analysis
