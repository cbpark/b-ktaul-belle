#ifndef B_KTAUL_BELLE_SRC_VARIABLE_H_
#define B_KTAUL_BELLE_SRC_VARIABLE_H_

#include <TRandom3.h>
#include <YAM2/yam2.h>
#include <memory>
#include <optional>
#include "input.h"
#include <utility>

namespace analysis {
class M2Reconstruction {
private:
    /// the M2 value (= -1 if the algorithm failed).
    double m2_;
    /// the reconstructed tau mass using the M2 solution (= -1 if no m2 solution
    /// found).
    double mtau_;

    M2Reconstruction(double m2, double mtau) : m2_(m2), mtau_(mtau) {}

public:
    M2Reconstruction() = delete;

    double m2() const { return m2_; }
    double mtau() const { return mtau_; }

    std::pair<double, double> get_result() const { return {m2_, mtau_}; }

    friend M2Reconstruction mkM2Reconstruction(
        const Input &input, const std::optional<yam2::M2Solution> &m2sol);
};

M2Reconstruction mkM2Reconstruction(
    const Input &input, const std::optional<yam2::M2Solution> &m2sol);

double mRecoilRandom(const Input &input, std::shared_ptr<TRandom3> rnd);
}  // namespace analysis

#endif  // B_KTAUL_BELLE_SRC_VARIABLE_H_
