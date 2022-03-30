#ifndef B_KTAUL_BELLE_SRC_INPUT_H_
#define B_KTAUL_BELLE_SRC_INPUT_H_

#include <Math/Vector3D.h>
#include <Math/Vector4D.h>

using LorentzVector = ROOT::Math::XYZTVector;
using Vector3 = ROOT::Math::XYZVectorF;

namespace analysis {
class Input {
private:
    /// Four-momentum of K_sig.
    LorentzVector k_sig_;
    /// Four-momentum of mu_sig.
    LorentzVector mu_sig_;

    Input(const LorentzVector &k_sig, const LorentzVector &mu_sig)
        : k_sig_(k_sig), mu_sig_(mu_sig) {}

public:
    Input() = delete;

    LorentzVector k_sig() const { return k_sig_; }
    LorentzVector mu_sig() const { return mu_sig_; }

    friend Input mkInputCM(const Vector3 &k_sig_v3, const Vector3 &mu_sig_v3);
};

Input mkInputCM(const Vector3 &k_sig_v3, const Vector3 &mu_sig_v3);
}  // namespace analysis

#endif  // B_KTAUL_BELLE_SRC_INPUT_H_
