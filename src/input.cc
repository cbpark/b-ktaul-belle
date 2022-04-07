#include "input.h"

#include <YAM2/yam2.h>
#include <cmath>  // std::sqrt
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
}  // namespace analysis
