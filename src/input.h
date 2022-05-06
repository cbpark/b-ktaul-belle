#ifndef B_KTAUL_BELLE_SRC_INPUT_H_
#define B_KTAUL_BELLE_SRC_INPUT_H_

#include <Math/Vector2D.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <YAM2/yam2.h>
#include <cmath>
#include <optional>
#include "constant.h"

using LorentzVector = ROOT::Math::XYZTVector;
using Vector3 = ROOT::Math::XYZVector;
using Vector3F = ROOT::Math::XYZVectorF;
using Vector2 = ROOT::Math::XYVector;
using Vector2F = ROOT::Math::XYVectorF;

namespace analysis {
class Input {
private:
    double sqrt_s_ = SQRTS;
    double pz_tot_ = PLONG;
    /// Four-momentum of K_sig.
    LorentzVector k_sig_;
    /// Four-momentum of mu_sig.
    LorentzVector mu_sig_;
    /// Four-momentum of htau_sig (hadronic tau decay products).
    LorentzVector htau_sig_;
    /// Four-momentum of D_tag.
    LorentzVector d_tag_;
    /// Four-Momentum of mu_tag.
    LorentzVector mu_tag_;

    /// missing transverse momentum.
    Vector2 ptmiss_;

    Vector3 vertex_bsig_;
    Vector3 vertex_btag_;

    // ptmiss is calculated by using visible particles.
    Input(const LorentzVector &k_sig, const LorentzVector &mu_sig,
          const LorentzVector &htau_sig, const LorentzVector &d_tag,
          const LorentzVector &mu_tag, const Vector3 &vertex_bsig,
          const Vector3 &vertex_btag)
        : k_sig_(k_sig),
          mu_sig_(mu_sig),
          htau_sig_(htau_sig),
          d_tag_(d_tag),
          mu_tag_(mu_tag),
          vertex_bsig_(vertex_bsig),
          vertex_btag_(vertex_btag) {
        LorentzVector pvis_tot = k_sig + mu_sig + htau_sig + d_tag + mu_tag;
        ptmiss_ = Vector2{-pvis_tot.px(), -pvis_tot.py()};
    }

    // ptmiss is an input.
    Input(const LorentzVector &k_sig, const LorentzVector &mu_sig,
          const LorentzVector &htau_sig, const LorentzVector &d_tag,
          const LorentzVector &mu_tag, const Vector2 &ptmiss,
          const Vector3 &vertex_bsig, const Vector3 &vertex_btag)
        : k_sig_(k_sig),
          mu_sig_(mu_sig),
          htau_sig_(htau_sig),
          d_tag_(d_tag),
          mu_tag_(mu_tag),
          ptmiss_(ptmiss),
          vertex_bsig_(vertex_bsig),
          vertex_btag_(vertex_btag) {}

public:
    Input() = delete;

    LorentzVector k_sig() const { return k_sig_; }
    LorentzVector mu_sig() const { return mu_sig_; }
    LorentzVector htau_sig() const { return htau_sig_; }
    LorentzVector d_tag() const { return d_tag_; }
    LorentzVector mu_tag() const { return mu_tag_; }

    double sqrt_s() const { return sqrt_s_; }
    void set_sqrt_s(double sqrt_s) { sqrt_s_ = sqrt_s; }

    double pz_tot() const { return pz_tot_; }
    void set_pz_tot(double pz_tot) { pz_tot_ = pz_tot; }

    // the energy of B meson at CM frame.
    double eb_cm() const { return 0.5 * sqrt_s_; }

    // |p| of B mesons at CM frame.
    double pb_cm() const { return std::sqrt(0.25 * sqrt_s_ * sqrt_s_ - MBSQ); }

    Vector2 ptmiss() const { return ptmiss_; }

    LorentzVector vis_sig() const { return k_sig_ + mu_sig_ + htau_sig_; }
    LorentzVector vis_tag() const { return d_tag_ + mu_tag_; }

    LorentzVector kl_sig() const { return k_sig_ + mu_sig_; }

    Vector3 vertex_bsig() const { return vertex_bsig_; }
    Vector3 vertex_btag() const { return vertex_btag_; }

    std::optional<yam2::InputKinematics> to_input_kinematics(
        double m_invisible) const;

    std::optional<yam2::InputKinematicsWithVertex>
    to_input_kinematics_with_vertex(
        std::optional<yam2::InputKinematics> &input_kinematics,
        double delta_theta) const;

    friend Input mkInput(const LorentzVector &k_sig,
                         const LorentzVector &mu_sig,
                         const LorentzVector &htau_sig,
                         const LorentzVector &d_tag,
                         const LorentzVector &mu_tag,
                         const std::optional<Vector2> &ptmiss,
                         const std::optional<Vector3> &ip_lab,
                         const std::optional<Vector3> &vertex_bsig_lab,
                         const std::optional<Vector3> &vertex_btag_lab);
};

Input mkInput(const LorentzVector &k_sig, const LorentzVector &mu_sig,
              const LorentzVector &htau_sig, const LorentzVector &d_tag,
              const LorentzVector &mu_tag,
              const std::optional<Vector2> &ptmiss = {},
              const std::optional<Vector3> &ip_lab = {},
              const std::optional<Vector3> &vertex_bsig_lab = {},
              const std::optional<Vector3> &vertex_btag_lab = {});
}  // namespace analysis

#endif  // B_KTAUL_BELLE_SRC_INPUT_H_
