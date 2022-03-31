#ifndef B_KTAUL_BELLE_SRC_INPUT_H_
#define B_KTAUL_BELLE_SRC_INPUT_H_

#include <Math/Vector2D.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>

using LorentzVector = ROOT::Math::XYZTVector;
using Vector3 = ROOT::Math::XYZVectorF;
using Vector2 = ROOT::Math::XYVector;

namespace analysis {
class Input {
private:
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

    Input(const LorentzVector &k_sig, const LorentzVector &mu_sig,
          const LorentzVector &htau_sig, const LorentzVector &d_tag,
          const LorentzVector &mu_tag)
        : k_sig_(k_sig),
          mu_sig_(mu_sig),
          htau_sig_(htau_sig),
          d_tag_(d_tag),
          mu_tag_(mu_tag) {
        LorentzVector pvis_tot = k_sig + mu_sig + htau_sig + d_tag + mu_tag;
        ptmiss_ = Vector2{-pvis_tot.px(), -pvis_tot.py()};
    }

public:
    Input() = delete;

    LorentzVector k_sig() const { return k_sig_; }
    LorentzVector mu_sig() const { return mu_sig_; }
    LorentzVector htau_sig() const { return htau_sig_; }
    LorentzVector d_tag() const { return d_tag_; }
    LorentzVector mu_tag() const { return mu_tag_; }
    Vector2 ptmiss() const { return ptmiss_; }

    friend Input mkInputCM(const Vector3 &k_sig_v3, const Vector3 &mu_sig_v3,
                           const Vector3 &htau_sig_v3, const Vector3 &d_tag_v3,
                           const Vector3 &mu_tag_v3);
};

Input mkInputCM(const Vector3 &k_sig_v3, const Vector3 &mu_sig_v3,
                const Vector3 &htau_sig_v3, const Vector3 &d_tag_v3,
                const Vector3 &mu_tag_v3);
}  // namespace analysis

#endif  // B_KTAUL_BELLE_SRC_INPUT_H_
