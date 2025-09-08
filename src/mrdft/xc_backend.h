#pragma once

#include <Eigen/Core>
#include <memory>
#include <string>

namespace mrdft {

class XCBackend {
public:
    virtual ~XCBackend() = default;

    // Configure functional & coefficient (same semantics as xcfun_set)
    virtual void set_functional(const std::string &name, double coeff) = 0;

    // Configure evaluation layout (order, spin, gamma vs explicit)
    virtual void configure(int order, bool spin, bool use_gamma) = 0;

    // Introspection
    virtual bool  is_gga() const = 0;
    virtual bool  is_metagga() const = 0;
    virtual double amount_exx() const = 0;
    virtual int   input_length() const = 0;
    virtual int   output_length() const = 0;

    // Evaluate at many grid points. Each ROW is one grid point.
    // Returns (nPts x nOut). Applies cutoff identical to old code.
    virtual Eigen::MatrixXd eval_transposed(const Eigen::MatrixXd &inp,
                                            double cutoff,
                                            bool is_spin_sep) const = 0;
};

using XCBackend_p = std::shared_ptr<XCBackend>;

// Stage 1 backend: wraps XCFun
XCBackend_p make_xcfun_backend();

} // namespace mrdft
