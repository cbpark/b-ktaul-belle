#ifndef B_KTAUL_BELLE_SRC_CONSTANT_H_
#define B_KTAUL_BELLE_SRC_CONSTANT_H_

#include <Math/Boost.h>
#include <Math/Vector4D.h>

namespace analysis {
constexpr double MK = 0.4937;

constexpr double MMU = 0.105;

constexpr double MPI = 0.13957;

constexpr double MD = 1.864;

constexpr double SQRTS = 10.583;

constexpr double EEM = 8.0;
constexpr double EEP = 3.5;
constexpr double EBEAMS = EEM + EEP;
constexpr double PLONG = EEM - EEP;

auto boostToCM() {
    ROOT::Math::XYZTVector sqrt_s{0.0, 0.0, PLONG, EBEAMS};
    return ROOT::Math::Boost(sqrt_s.BoostToCM());
}
}  // namespace analysis

#endif  // B_KTAUL_BELLE_SRC_CONSTANT_H_
