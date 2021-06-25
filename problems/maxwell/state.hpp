#ifndef PROBLEMS_MAXWELL_STATE_HPP_
#define PROBLEMS_MAXWELL_STATE_HPP_

#include "ads/lin/tensor.hpp"


namespace ads {

struct state {
    using field = ads::lin::tensor<double, 3>;

    field E1;
    field E2;
    field E3;

    field H1;
    field H2;
    field H3;

    state(const std::array<std::size_t, 3> shape)
    : E1{shape}
    , E2{shape}
    , E3{shape}
    , H1{shape}
    , H2{shape}
    , H3{shape}
    { }

    void clear() {
        zero(E1);
        zero(E2);
        zero(E3);

        zero(H1);
        zero(H2);
        zero(H3);
    }
};

}

#endif // ifndef PROBLEMS_MAXWELL_STATE_HPP_