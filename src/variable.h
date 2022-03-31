#ifndef B_KTAUL_BELLE_SRC_VARIABLE_H_
#define B_KTAUL_BELLE_SRC_VARIABLE_H_

#include <TRandom3.h>
#include <memory>
#include "input.h"

namespace analysis {
double mRecoilRandom(const Input &input, std::shared_ptr<TRandom3> rnd);
}  // namespace analysis

#endif  // B_KTAUL_BELLE_SRC_VARIABLE_H_
