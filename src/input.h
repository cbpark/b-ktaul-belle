#ifndef B_KTAUL_BELLE_SRC_INPUT_H_
#define B_KTAUL_BELLE_SRC_INPUT_H_

#include <Math/Vector2D.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <YAM2/yam2.h>
#include <optional>
#include "constant.h"

using LorentzVector = ROOT::Math::XYZTVector;
using Vector3F = ROOT::Math::XYZVectorF;
using Vector2 = ROOT::Math::XYVector;
using Vector2F = ROOT::Math::XYVectorF;

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

    Input(const LorentzVector &k_sig, const LorentzVector &mu_sig,
          const LorentzVector &htau_sig, const LorentzVector &d_tag,
          const LorentzVector &mu_tag, const Vector2 &ptmiss)
        : k_sig_(k_sig),
          mu_sig_(mu_sig),
          htau_sig_(htau_sig),
          d_tag_(d_tag),
          mu_tag_(mu_tag),
          ptmiss_(ptmiss) {}

public:
    Input() = delete;

    LorentzVector k_sig() const { return k_sig_; }
    LorentzVector mu_sig() const { return mu_sig_; }
    LorentzVector htau_sig() const { return htau_sig_; }
    LorentzVector d_tag() const { return d_tag_; }
    LorentzVector mu_tag() const { return mu_tag_; }
    Vector2 ptmiss() const { return ptmiss_; }

    LorentzVector vis_sig() const { return k_sig_ + mu_sig_ + htau_sig_; }
    LorentzVector vis_tag() const { return d_tag_ + mu_tag_; }

    LorentzVector kl_sig() const { return k_sig_ + mu_sig_; }

    std::optional<yam2::InputKinematics> to_input_kinematics(
        double m_invisible, double pz_tot) const;

    friend Input mkInputCM(const LorentzVector &k_sig,
                           const LorentzVector &mu_sig,
                           const LorentzVector &htau_sig,
                           const LorentzVector &d_tag,
                           const LorentzVector &mu_tag,
                           const std::optional<Vector2> &ptmiss);

    // template <typename V3>
    // friend Input mkInputCM(const V3 &k_sig_v3, const V3 &mu_sig_v3,
    //                        const V3 &htau_sig_v3, const V3 &d_tag_v3,
    //                        const V3 &mu_tag_v3);

    // template <typename V3>
    // friend Input mkInputCM(const Vector2F &b_sig_v2, const Vector2F
    // &b_tag_v2,
    //                        const V3 &k_sig_v3, const V3 &mu_sig_v3,
    //                        const V3 &htau_sig_v3, const V3 &d_tag_v3,
    //                        const V3 &mu_tag_v3);
};

// template <typename V3>
// LorentzVector toLorentzVector(const V3 &v3, double mass) {
//     double energy = std::sqrt(v3.mag2() + mass * mass);
//     return {v3.x(), v3.y(), v3.z(), energy};
// }

inline Input mkInputCM(const LorentzVector &k_sig, const LorentzVector &mu_sig,
                       const LorentzVector &htau_sig,
                       const LorentzVector &d_tag, const LorentzVector &mu_tag,
                       const std::optional<Vector2> &ptmiss = {}) {
    // std::cout << "k_sig(before): " << k_sig << '\n';
    LorentzVector k_sig_cm = boostToCM()(k_sig);
    // std::cout << "mu_sig(before): " << mu_sig << '\n';
    LorentzVector mu_sig_cm = boostToCM()(mu_sig);
    // std::cout << "htau_sig(before): " << htau_sig << '\n';
    LorentzVector htau_sig_cm = boostToCM()(htau_sig);
    // std::cout << "d_tag(before): " << d_tag << '\n';
    LorentzVector d_tag_cm = boostToCM()(d_tag);
    // std::cout << "mu_tag(before): " << mu_tag << '\n';
    LorentzVector mu_tag_cm = boostToCM()(mu_tag);

    if (!ptmiss) {
        return {k_sig_cm, mu_sig_cm, htau_sig_cm, d_tag_cm, mu_tag_cm};
    } else {
        return {k_sig_cm, mu_sig_cm, htau_sig_cm,
                d_tag_cm, mu_tag_cm, ptmiss.value()};
    }
}

// template <typename V3>
// Input mkInputCM(const V3 &k_sig_v3, const V3 &mu_sig_v3, const V3
// &htau_sig_v3,
//                 const V3 &d_tag_v3, const V3 &mu_tag_v3) {
//     LorentzVector k_sig = toLorentzVector(k_sig_v3, MK);
//     // std::cout << "k_sig(before): " << k_sig << '\n';
//     k_sig = boostToCM()(k_sig);
//     // std::cout << "k_sig(after): " << k_sig << '\n';

//     LorentzVector mu_sig = toLorentzVector(mu_sig_v3, MMU);
//     // std::cout << "mu_sig(before): " << mu_sig << '\n';
//     mu_sig = boostToCM()(mu_sig);
//     // std::cout << "mu_sig(after): " << mu_sig << '\n';

//     LorentzVector htau_sig = toLorentzVector(htau_sig_v3, MPI);
//     // std::cout << "htau_sig(before): " << htau_sig << '\n';
//     htau_sig = boostToCM()(htau_sig);
//     // std::cout << "htau_sig(after): " << htau_sig << '\n';

//     LorentzVector d_tag = toLorentzVector(d_tag_v3, MD);
//     // std::cout << "d_tag(before): " << d_tag << '\n';
//     d_tag = boostToCM()(d_tag);
//     // std::cout << "d_tag(after): " << d_tag << '\n';

//     LorentzVector mu_tag = toLorentzVector(mu_tag_v3, MMU);
//     // std::cout << "mu_tag(before): " << mu_tag << '\n';
//     mu_tag = boostToCM()(mu_tag);
//     // std::cout << "mu_tag(after): " << mu_tag << '\n';

//     return {k_sig, mu_sig, htau_sig, d_tag, mu_tag};
// }

// template <typename V3>
// Input mkInputCM(const Vector2F &b_sig_v2, const Vector2F &b_tag_v2,
//                 const V3 &k_sig_v3, const V3 &mu_sig_v3, const V3
//                 &htau_sig_v3, const V3 &d_tag_v3, const V3 &mu_tag_v3) {
//     LorentzVector k_sig = toLorentzVector(k_sig_v3, MK);
//     // std::cout << "k_sig(before): " << k_sig << '\n';
//     k_sig = boostToCM()(k_sig);
//     // std::cout << "k_sig(after): " << k_sig << '\n';

//     LorentzVector mu_sig = toLorentzVector(mu_sig_v3, MMU);
//     // std::cout << "mu_sig(before): " << mu_sig << '\n';
//     mu_sig = boostToCM()(mu_sig);
//     // std::cout << "mu_sig(after): " << mu_sig << '\n';

//     LorentzVector htau_sig = toLorentzVector(htau_sig_v3, MPI);
//     // std::cout << "htau_sig(before): " << htau_sig << '\n';
//     htau_sig = boostToCM()(htau_sig);
//     // std::cout << "htau_sig(after): " << htau_sig << '\n';

//     LorentzVector d_tag = toLorentzVector(d_tag_v3, MD);
//     // std::cout << "d_tag(before): " << d_tag << '\n';
//     d_tag = boostToCM()(d_tag);
//     // std::cout << "d_tag(after): " << d_tag << '\n';

//     LorentzVector mu_tag = toLorentzVector(mu_tag_v3, MMU);
//     // std::cout << "mu_tag(before): " << mu_tag << '\n';
//     mu_tag = boostToCM()(mu_tag);
//     // std::cout << "mu_tag(after): " << mu_tag << '\n';

//     // We don't need to boost the B mesons because the transverse components
//     do
//     // not transform.
//     auto pt_bb = b_sig_v2 + b_tag_v2;
//     LorentzVector pvis_tot = k_sig + mu_sig + htau_sig + d_tag + mu_tag;
//     Vector2 ptmiss{pt_bb.x() - pvis_tot.px(), pt_bb.y() - pvis_tot.py()};

//     return {k_sig, mu_sig, htau_sig, d_tag, mu_tag, ptmiss};
// }
}  // namespace analysis

#endif  // B_KTAUL_BELLE_SRC_INPUT_H_
