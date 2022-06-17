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

yam2::SpatialMomentum fromVector3(const Vector3 &p3) {
    return {p3.x(), p3.y(), p3.z()};
}

std::optional<yam2::InputKinematics> Input::to_input_kinematics(
    double m_invisible1, double m_invisible2) const {
    auto p1 = fromLorentzVector(this->vis_sig());
    auto p2 = fromLorentzVector(this->vis_tag());

    auto ptmiss = fromVector2(this->ptmiss());
    yam2::Mass m_parent{MB};

    return yam2::mkInput(p1, p2, ptmiss, yam2::Mass{m_invisible1},
                         yam2::Mass{m_invisible2}, m_parent, m_parent, sqrt_s_,
                         {pz_tot_});
}

std::optional<yam2::InputKinematicsWithVertex>
Input::to_input_kinematics_with_vertex(
    std::optional<yam2::InputKinematics> &input_kinematics,
    double delta_theta) const {
    auto vertex1 = fromVector3(this->vertex_bsig());
    auto vertex2 = fromVector3(this->vertex_btag());

    return yam2::mkInputWithVertex(input_kinematics, vertex1, vertex2,
                                   delta_theta);
}

Input mkInput(const LorentzVector &k_sig, const LorentzVector &mu_sig,
              const LorentzVector &htau_sig, const LorentzVector &d_tag,
              const LorentzVector &mu_tag, const std::optional<Vector2> &ptmiss,
              const std::optional<Vector3> &ip_lab,
              const std::optional<Vector3> &vertex_bsig_lab,
              const std::optional<Vector3> &vertex_btag_lab) {
    Vector3 vertex_bsig{0.0, 0.0, 0.0};
    Vector3 vertex_btag{0.0, 0.0, 0.0};
    if (ip_lab && vertex_bsig_lab && vertex_btag_lab) {
        vertex_bsig = vertex_bsig_lab.value() - ip_lab.value();
        vertex_btag = vertex_btag_lab.value() - ip_lab.value();
    }

    if (!ptmiss) {
        return {k_sig,  mu_sig,      htau_sig,   d_tag,
                mu_tag, vertex_bsig, vertex_btag};
    } else {
        return {k_sig,  mu_sig,         htau_sig,    d_tag,
                mu_tag, ptmiss.value(), vertex_bsig, vertex_btag};
    }
}
}  // namespace analysis
