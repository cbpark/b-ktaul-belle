#include "variable.h"

#include <TRandom3.h>
#include <cmath>  // std::sqrt
#include <memory>
#include "constant.h"
#include "input.h"

namespace analysis {
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

double mTauM2s(const Input &input) {
    return 0.0;
}
}  // namespace analysis
