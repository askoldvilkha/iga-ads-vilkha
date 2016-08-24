#ifndef PROBLEMS_ELASTICITY_ELASTICITY_HPP_
#define PROBLEMS_ELASTICITY_ELASTICITY_HPP_

#include <cmath>
#include "ads/simulation.hpp"
#include "ads/executor/galois.hpp"
#include "ads/output_manager.hpp"


namespace problems {

    class linear_elasticity : public ads::simulation_3d {
        using Base = simulation_3d;
        using tensor = double[3][3];

        struct state {
            vector_type ux, uy, uz;
            vector_type vx, vy, vz;
            state(std::array<std::size_t, 3> shape)
            : ux{ shape }, uy{ shape }, uz{ shape }
            , vx{ shape }, vy{ shape }, vz{ shape }
            { }
        };

        state now, prev;

        ads::output_manager<3> output;

        ads::galois_executor executor{8};

        template <typename Fun>
        void for_all(state& s, Fun fun) {
            fun(s.ux);
            fun(s.uy);
            fun(s.uz);
            fun(s.vx);
            fun(s.vx);
            fun(s.vx);
        }

    public:
        linear_elasticity(const ads::config_3d& config)
        : Base{ config }
        , now{ shape() }, prev{ shape() }
        , output{ x.B, y.B, z.B, 50 }
        { }

    private:

        void before() override {
            prepare_matrices();

            // auto init = [this](double x, double y, double z) { return x; };
            // projection(now.ux, init);
            // solve(now.ux);

            output.to_file("init.vti",
                           output.evaluate(now.ux),
                           output.evaluate(now.uy),
                           output.evaluate(now.uz));
        }

        void compute_rhs(double t) {
            for_all(now, [](vector_type& a) { zero(a); });
            executor.for_each(elements(), [&](index_type e) {
                auto local = local_contribution(e, t);
                executor.synchronized([&] {
                    apply_local_contribution(local, e);
                });
            });
        }

        state local_contribution(index_type e, double t) const {
            auto local = state{ local_shape() };
            double J = jacobian(e);
            for (auto q : quad_points()) {
                auto x = point(e, q);
                double w = weigth(q);
                for (auto a : dofs_on_element(e)) {
                    value_type ux = eval_fun(prev.ux, e, q);
                    value_type uy = eval_fun(prev.uy, e, q);
                    value_type uz = eval_fun(prev.uz, e, q);
                    value_type vx = eval_fun(prev.vx, e, q);
                    value_type vy = eval_fun(prev.vy, e, q);
                    value_type vz = eval_fun(prev.vz, e, q);

                    tensor eps = {
                        {         ux.dx,         0.5 * (ux.dy + uy.dx), 0.5 * (ux.dz + uz.dx) },
                        { 0.5 * (ux.dy + uy.dx),         uy.dy,         0.5 * (uy.dz + uz.dy) },
                        { 0.5 * (ux.dz + uz.dx), 0.5 * (uy.dz + uz.dy),         uz.dz         }
                    };
                    tensor s;
                    stress_tensor(s, eps);

                    value_type b = eval_basis(e, q, a);

                    auto F = force(x, t);
                    double rho = 1;
                    double axb = (-eps[0][0] * b.dx - eps[0][1] * b.dy - eps[0][2] * b.dz + F[0] * b.val) / rho;
                    double ayb = (-eps[1][0] * b.dx - eps[1][1] * b.dy - eps[1][2] * b.dz + F[1] * b.val) / rho;
                    double azb = (-eps[2][0] * b.dx - eps[2][1] * b.dy - eps[2][2] * b.dz + F[2] * b.val) / rho;

                    double t2 = steps.dt * steps.dt / 2;
                    auto aa = dof_global_to_local(e, a);
                    ref(local.ux, aa) += ((ux.val + steps.dt * vx.val) * b.val + t2 * axb) * w * J;
                    ref(local.uy, aa) += ((uy.val + steps.dt * vy.val) * b.val + t2 * ayb) * w * J;
                    ref(local.uz, aa) += ((uz.val + steps.dt * vz.val) * b.val + t2 * azb) * w * J;

                    ref(local.vx, aa) += (ux.val * b.val + steps.dt * axb) * w * J;
                    ref(local.vy, aa) += (uy.val * b.val + steps.dt * ayb) * w * J;
                    ref(local.vz, aa) += (uz.val * b.val + steps.dt * azb) * w * J;
                }
            }
            return local;
        }

        void apply_local_contribution(const state& loc, index_type e) {
            update_global_rhs(now.ux, loc.ux, e);
            update_global_rhs(now.uy, loc.uy, e);
            update_global_rhs(now.uz, loc.uz, e);
            update_global_rhs(now.vx, loc.vx, e);
            update_global_rhs(now.vy, loc.vy, e);
            update_global_rhs(now.vz, loc.vz, e);
        }


        void compute_rhs__(double t) {
            for (auto e : elements()) {
                double J = jacobian(e);
                for (auto q : quad_points()) {
                    auto x = point(e, q);
                    double w = weigth(q);
                    for (auto a : dofs_on_element(e)) {
                        value_type ux = eval_fun(prev.ux, e, q);
                        value_type uy = eval_fun(prev.uy, e, q);
                        value_type uz = eval_fun(prev.uz, e, q);
                        value_type vx = eval_fun(prev.vx, e, q);
                        value_type vy = eval_fun(prev.vy, e, q);
                        value_type vz = eval_fun(prev.vz, e, q);

                        tensor eps = {
                            {         ux.dx,         0.5 * (ux.dy + uy.dx), 0.5 * (ux.dz + uz.dx) },
                            { 0.5 * (ux.dy + uy.dx),         uy.dy,         0.5 * (uy.dz + uz.dy) },
                            { 0.5 * (ux.dz + uz.dx), 0.5 * (uy.dz + uz.dy),         uz.dz         }
                        };
                        tensor s{};
                        stress_tensor(s, eps);

                        value_type b = eval_basis(e, q, a);

                        auto F = force(x, t);
                        double rho = 1;
                        double axb = (-eps[0][0] * b.dx - eps[0][1] * b.dy - eps[0][2] * b.dz + F[0] * b.val) / rho;
                        double ayb = (-eps[1][0] * b.dx - eps[1][1] * b.dy - eps[1][2] * b.dz + F[1] * b.val) / rho;
                        double azb = (-eps[2][0] * b.dx - eps[2][1] * b.dy - eps[2][2] * b.dz + F[2] * b.val) / rho;

                        double t2 = steps.dt * steps.dt / 2;
                        ref(now.ux, a) += ((ux.val + steps.dt * vx.val) * b.val + t2 * axb) * w * J;
                        ref(now.uy, a) += ((uy.val + steps.dt * vy.val) * b.val + t2 * ayb) * w * J;
                        ref(now.uz, a) += ((uz.val + steps.dt * vz.val) * b.val + t2 * azb) * w * J;

                        ref(now.vx, a) += (ux.val * b.val + steps.dt * axb) * w * J;
                        ref(now.vy, a) += (uy.val * b.val + steps.dt * ayb) * w * J;
                        ref(now.vz, a) += (uz.val * b.val + steps.dt * azb) * w * J;
                    }
                }
            }
        }

        void stress_tensor(tensor& s, const tensor& eps) const {
            double lambda = 0.1;
            double mi = 0.1;
            double tr = eps[0][0] + eps[1][1] + eps[2][2];

            for (int i = 0; i < 3; ++ i) {
                for (int j = 0; j < 3; ++ j) {
                    s[i][j] = 2 * mi * eps[i][j];
                }
                s[i][i] += lambda * tr;
            }
        }

        point_type force(point_type x, double t) const {
            using std::pow;
            constexpr double t0 = 0.02;
            double tt = t / t0;
            double f = pow(tt * (1 - tt), 2);
            double r = pow(x[0] - 1, 2) + pow(x[1] - 1, 2) + pow(x[2] - 1, 2);
            double a = - 10 * f * std::exp(- 10 * r);
            return {a, a, a};
        }

        void before_step(int /*iter*/, double /*t*/) override {
            using std::swap;
            swap(now, prev);
        }

        void step(int /*iter*/, double t) override {
            compute_rhs(t);
            for_all(now, [this](vector_type& a) { solve(a); });
        }

        void after_step(int iter, double /*t*/) override {
            std::cout << "Iteration " << iter << std::endl;
            if (iter % 10 == 0) {
                std::cout << "Outputting..." << std::endl;
                output.to_file("out_%d.vti", iter,
                               output.evaluate(now.ux),
                               output.evaluate(now.uy),
                               output.evaluate(now.uz));
            }
        }


        double& ref(vector_type& v, index_type idx) const {
            return v(idx[0], idx[1], idx[2]);
        }
    };
}



#endif /* PROBLEMS_ELASTICITY_ELASTICITY_HPP_ */