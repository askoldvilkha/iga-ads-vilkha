// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef CO2_SEQUESTRATION_CO2_SEQUESTRATION_2D_HPP
#define CO2_SEQUESTRATION_CO2_SEQUESTRATION_2D_HPP

#include <galois/Timer.h>

#include "ads/executor/galois.hpp"
#include "ads/output_manager.hpp"
#include "ads/simulation.hpp"
#include "ads/solver/mumps.hpp"

#include "ads/lin/tensor/tensor.hpp"

namespace ads::problems {

class co2_sequestration_2d : public simulation_2d {
private:
    using Base = simulation_2d;
    vector_type p, s, p_prev, s_prev;

    output_manager<2> output;
    galois_executor executor;
    galois::StatTimer initialization_timer{"initialization"};
    galois::StatTimer integration_timer_s{"integration_s"};
    galois::StatTimer integration_timer_p{"integration_p"};
    galois::StatTimer solver_timer_p{"solver_p"};
    galois::StatTimer solver_timer_s{"solver_s"};
    galois::StatTimer output_timer{"output"};
    ads::mumps::solver solver;

    std::ofstream normL2_file{"L2_norms.dat"};

    // parameter maps
    using map_array = ads::lin::tensor<double, 2>;
    map_array porosity;
    map_array permeability;

    // parameters
    const double mu_w;     // brine viscosity
    const double mu_g;     // gas viscosity
    const double K;        // permeability tensor
    const double phi;      // porosity
    const double rho_w;    // brine density
    const double rho_g;    // gas density
    const double g;        // gravitational acceleration
    const double qg_x;     // gas injection location - x
    const double qg_y;     // gas injection location - y
    const double qg_rate;  // gas injection rate
    const double qg_radius; // gas injection radius (squared)
    const int mesh_x;
    const int mesh_y;
    std::string porosity_map;
    std::string permeability_map;
    const double darcy_to_si = 9.869233e-13; // conversion factor from Darcy to m^2
    bool verbose;

public:
    explicit co2_sequestration_2d(const config_2d& config, const double mu_w, const double mu_g,
                                  const double K, const double phi,
                                  const double rho_w, const double rho_g, const double g,
                                  const double qg_x, const double qg_y, const double qg_rate, const double qg_radius,
                                  std::string porosity_map, std::string permeability_map, bool verbose, int threads)
    : Base{config}
    , p{shape()}
    , s{shape()}
    , p_prev{shape()}
    , s_prev{shape()}
    , mu_w{mu_w}
    , mu_g{mu_g}
    , K{K}
    , phi{phi}
    , rho_w{rho_w}
    , rho_g{rho_g}
    , g{g}
    , qg_x{qg_x}
    , qg_y{qg_y}
    , qg_rate{qg_rate}
    , qg_radius{qg_radius}
    , verbose{verbose}
    , porosity_map{porosity_map}
    , permeability_map{permeability_map}
    , mesh_x{config.x.b}
    , mesh_y{config.y.b}
    , output{x.B, y.B, config.x.elements, config.y.elements}
    , porosity{{1, 1}}
    , permeability{{1, 1}}
    , executor{threads} { }

    // this function sets the initial state of the gas saturation
    double init_state(double x, double y) {
        /*double dx = x - 50;
        double dy = y - 8;
        double r2 = std::min(0.25 * (dx * dx + dy * dy), 1.0);
        return 0.5 * ((r2 - 1) * (r2 - 1) * (r2 + 1) * (r2 + 1));
        */
        return 0;
    };

    double source_g(double x, double y, double t) {
        double dx = x - qg_x;
        double dy = y - qg_y;
        double r2 = std::min(1.0 / qg_radius * (dx * dx + dy * dy), 1.0);
        return qg_rate * ((r2 - 1) * (r2 - 1) * (r2 + 1) * (r2 + 1));
    }

    double source_w(double x, double y, double t) {
        double val = 0;
        return val;
    }

private:
    auto read_map_data(std::string_view filename) -> ads::lin::tensor<double, 2> {
        std::ifstream input{std::string(filename)};

        int nx;
        int ny;
        input >> nx >> ny;
        auto data = ads::lin::tensor<double, 2>{{nx, ny}};

        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                double x;
                double y;
                double val;
                input >> x >> y >> val;
                data(i, j) = val;
            }
        }
        return data;
    }

    static auto point_lookup(double x, int n) -> std::tuple<int, int> {
        int min_1 = std::min(static_cast<int>(x * n), n - 1);
        int min_2 = std::min(static_cast<int>(min_1 + 1), n);

        if (min_1 == min_2) {
            min_2 = min_1 + 1;
            std::cout << "min_1 == min_2" << std::endl;
        }

        return std::make_tuple(min_1, min_2);
    }

    auto approximate_map_data(point_type xy, map_array map) -> double {
        double x = xy[0];
        double y = xy[1];
        double mesh_xf = mesh_x;
        double mesh_yf = mesh_y;
        const auto [ix1, ix2] = point_lookup(x / mesh_xf, map.size(0) - 1);
        const auto [iy1, iy2] = point_lookup(y / mesh_yf, map.size(1) - 1);

        // we will use bilinear interpolation here
        double f_x1y1 = map(ix1, iy1);
        double f_x2y1 = map(ix2, iy1);
        double f_x1y2 = map(ix1, iy2);
        double f_x2y2 = map(ix2, iy2);

        double x1 = ix1 * mesh_xf / (map.size(0) - 1);
        double x2 = ix2 * mesh_xf / (map.size(0) - 1);
        double y1 = iy1 * mesh_yf / (map.size(1) - 1);
        double y2 = iy2 * mesh_yf / (map.size(1) - 1);

        double fx_y1 = f_x1y1 + (f_x2y1 - f_x1y1) * (x - x1) / (x2 - x1);
        double fx_y2 = f_x1y2 + (f_x2y2 - f_x1y2) * (x - x1) / (x2 - x1);
        double fxy = fx_y1 + (fx_y2 - fx_y1) * (y - y1) / (y2 - y1);

        if (fxy < 1e-4 || std::isnan(fxy)) {
            std::cout << x << " " << y << std::endl;
            std::cout << ix1 << " " << ix2 << " " << iy1 << " " << iy2 << std::endl;
            std::cout << "Interpolation function testing at: x = " << x << ", y = " << y << std::endl;
            std::cout << "x1: " << x1 << ", x2: " << x2 << ", y1: " << y1 << ", y2: " << y2 << std::endl;
            std::cout << "f_x1y1: " << f_x1y1 << ", f_x2y1: " << f_x2y1 << ", f_x1y2: " << f_x1y2
                      << ", f_x2y2: " << f_x2y2 << std::endl;
            std::cout << "fx_y1: " << fx_y1 << ", fx_y2: " << fx_y2 << std::endl;
            std::cout << "fxy: " << fxy << std::endl;
            std::cout << "map size: " << map.size(0) << ", " << map.size(1) << std::endl;
            std::cout << std::endl;
        }

        return fxy;
    }

    void visualize_map_data(map_array map, const std::string& output_file) {

        std::ofstream ofs(output_file);
        if (ofs.is_open()) {
            for (int i = 0; i < mesh_x; i++) {
                for (int j = 0; j < mesh_y * 4; j++) {
                    double x = i * mesh_x / (mesh_x);
                    double y = j * mesh_y / (mesh_y * 4);
                    double val = approximate_map_data({x, y}, map);
                    ofs << x << " " << y << " " << val << std::endl;
                }
            }
        }
    }

    void solve(vector_type& v) {
        // lin::vector buf{{y.dofs()}};
        // compute_projection(buf, y.basis, [](double y) { return std::sin(y * M_PI); });
        // for (int i = 0; i < y.dofs(); ++i) {
        //    v(0, i) = buf(i);
        // }
        Base::solve(v);
    }

    void prepare_matrices() {
        x.fix_left();
        x.fix_right();
        y.fix_left();
        y.fix_right();

        // !!all the bcs are Dirichlet for now - this is not correct for the final version!!
        Base::prepare_matrices();
    }

    void before() override {
        initialization_timer.start();
        prepare_matrices();

        if (porosity_map != "none") {
            porosity = read_map_data(porosity_map);
            visualize_map_data(porosity, "porosity_visualization.dat");
        }

        if (permeability_map != "none") {
            permeability = read_map_data(permeability_map);
            visualize_map_data(permeability, "permeability_visualization.dat");
        }

        auto init = [this](double x, double y) { return init_state(x, y); };

        projection(s, init);
        solve(s);

        output.to_file(p, "p.init.data");
        output.to_file(s, "s.init.data");
        if (verbose) {
            std::cout << "Initial projection computed" << std::endl;
        }
        initialization_timer.stop();
    }

    void before_step(int /*iter*/, double /*t*/) override {
        using std::swap;
        swap(p, p_prev);
        swap(s, s_prev);
    }

    void step(int /*iter*/, double t) override {
        // solve for p
        integration_timer_p.start();
        compute_rhs_p(t);
        // compute_rhs_simple(t);
        // p(0, 0) = 0;

        dirichlet_bc(p, boundary::left, x, y, [](double t) { return 0; });
        dirichlet_bc(p, boundary::right, x, y, [](double t) { return 0; });
        dirichlet_bc(p, boundary::bottom, x, y, [](double t) { return 0; });
        dirichlet_bc(p, boundary::top, x, y, [](double t) { return 0; });

        ads::mumps::problem problem_p(p.data(), p.size());
        assemble_problem(problem_p);
        integration_timer_p.stop();

        solver_timer_p.start();
        solver.solve(problem_p);
        solver_timer_p.stop();

        // once p is solved, we can solve for s and move to the next iteration afterwards

        integration_timer_s.start();
        compute_rhs(t);
        dirichlet_bc(s, boundary::left, x, y, [](double t) { return 0; });
        dirichlet_bc(s, boundary::right, x, y, [](double t) { return 0; });
        dirichlet_bc(s, boundary::bottom, x, y, [](double t) { return 0; });
        dirichlet_bc(s, boundary::top, x, y, [](double t) { return 0; });
        integration_timer_s.stop();

        solver_timer_s.start();
        solve(s);
        solver_timer_s.stop();
    }

    void after_step(int iter, double /*t*/) override {
        if (iter % 10 == 0) {
            output_timer.start();
            output.to_file(p, "p.out_%d.data", iter);
            output.to_file(s, "s.out_%d.data", iter);
            normL2_file << iter << " " << normL2(s, x, y) <<  " " << normL2(p, x, y) << std::endl;
            output_timer.stop();
            if (verbose) {
                std::cout << "Iteration " << iter << " passed" << std::endl;
            }
        }
    }

    void after() override {
        std::cout << "\n=== Timer Results ===" << std::endl;
        std::cout << "Initialization: " << initialization_timer.get() << " ms" << std::endl;
        std::cout << "Integration (s): " << integration_timer_s.get() << " ms" << std::endl;
        std::cout << "Integration (p): " << integration_timer_p.get() << " ms" << std::endl;
        std::cout << "Solver (p): " << solver_timer_p.get() << " ms" << std::endl;
        std::cout << "Solver (s): " << solver_timer_s.get() << " ms" << std::endl;
        std::cout << "Output: " << output_timer.get() << " ms" << std::endl;
        std::cout << "=====================\n" << std::endl;
        std::cout << "Final L2 norms: " << std::endl;
        std::cout << "||p||_L2 = " << normL2(p, x, y) << std::endl;
        std::cout << "||s||_L2 = " << normL2(s, x, y) << std::endl;
    }

    void compute_rhs_simple(double t) {
        auto& rhs = p;

        zero(rhs);

        executor.for_each(elements(), [&](index_type e) {
            auto U = element_rhs();

            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weight(q);
                auto x = point(e, q);
                for (auto a : dofs_on_element(e)) {
                    auto aa = dof_global_to_local(e, a);
                    value_type v = eval_basis(e, q, a);

                    // double const n = 8.0 / 100.0;
                    // double const m = 11.0 / 16.0;
                    // auto const pi = M_PI;

                    // auto const lambda = pi * pi * (n * n + m * m);
                    auto fval = source_g(x[0], x[1], t);

                    double val = fval * v.val;
                    U(aa[0], aa[1]) += val * w * J;
                }
            }

            executor.synchronized([&]() { update_global_rhs(rhs, U, e); });
        });
    }

    double approximate_K_at_point(point_type xy) {
        double K_here;
        if (permeability_map != "none") {
            K_here = approximate_map_data(xy, permeability);
            K_here = K_here * darcy_to_si * 1e-3; // convert from mD to m^2
        }
        else {
            K_here = K;
        }
        return K_here;
    }

    void compute_rhs(double t) {
        auto& rhs = s;

        zero(rhs);

        executor.for_each(elements(), [&](index_type e) {
            auto U = element_rhs();

            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weight(q);
                auto x = point(e, q);
                for (auto a : dofs_on_element(e)) {
                    auto aa = dof_global_to_local(e, a);
                    value_type v = eval_basis(e, q, a);
                    value_type s = eval_fun(s_prev, e, q);
                    value_type p_here = eval_fun(p, e, q);

                    double K_here = approximate_K_at_point(x);
                    double s_val = std::clamp(s.val, 0.0, 1.0);
                    double term_1 = s_val * grad_dot(p_here, v) * K_here / mu_g;
                    double term_2 = s_val * v.dy * K_here * rho_g * g / mu_g;
                    double term_3 = v.val * source_g(x[0], x[1], t);

                    double phi_here;
                    if (porosity_map != "none") {
                        phi_here = approximate_map_data(x, porosity);
                    } else {
                        phi_here = phi;
                    }

                    // temporary porosity increase
                    phi_here += 0.1;
                    double val = (term_2 - term_1 + term_3) * steps.dt / phi_here + s_val * v.val;

                    // NOTE! this term is a temporary enforcement of the upper ceiling on saturation
                    double term_extra = -1 * s.val * v.val * (x[1] >= mesh_y - 1);
                    val += term_extra;

                    U(aa[0], aa[1]) += val * w * J;
                }
            }

            executor.synchronized([&]() { update_global_rhs(rhs, U, e); });
        });
    }

    void compute_rhs_p(double t) {
        auto& rhs = p;

        zero(rhs);

        executor.for_each(elements(), [&](index_type e) {
            auto U = element_rhs();

            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weight(q);
                auto x = point(e, q);
                for (auto a : dofs_on_element(e)) {
                    auto aa = dof_global_to_local(e, a);
                    value_type v = eval_basis(e, q, a);
                    value_type s = eval_fun(s_prev, e, q);

                    double K_here = approximate_K_at_point(x);
                    //double term_1 = K_here * s.val * v.dy * g * rho_w / mu_w;
                    // considering no differentiation in weak form
                    double term_1 = K_here * (1 - s.val) * v.dy * g * rho_w / mu_w;
                    double term_2 = K_here * s.val * v.dy * g * rho_g / mu_g;
                    double term_3 = (source_w(x[0], x[1], t) + source_g(x[0], x[1], t)) * v.val;

                    // double val = term_1 - term_2 - term_3;
                    // considering negative sign in front of divergence in weak form
                    double val = term_1 + term_2 + term_3;
                    U(aa[0], aa[1]) += val * w * J;
                }
            }

            executor.synchronized([&]() { update_global_rhs(rhs, U, e); });
        });
    }

    void assemble_problem(ads::mumps::problem& problem) {
        executor.for_each(dofs(x, y), [&](auto a) {
            std::vector<std::tuple<int, int, double>> vals_buf;
            for (auto b : overlapping_dofs(a, x, y)) {
                if (is_fixed(a, x, y))
                    continue;

                double val = 0;
                for (auto e : elements_supporting_dof(a, x, y)) {
                    if (!supported_in(b, e, x, y))
                        continue;

                    double J = jacobian(e, x, y);
                    for (auto q : quad_points(x, y)) {
                        double w = weight(q, x, y);
                        value_type ww = eval_basis(e, q, a, x, y);
                        value_type uu = eval_basis(e, q, b, x, y);

                        double s = std::clamp(eval_fun(s_prev, e, q).val, 0.0, 1.0);
                        double diff_1 = (1 - s) / mu_w;
                        double diff_2 = s / mu_g;

                        double K_here = approximate_K_at_point(point(e, q, x, y));

                        // double bwu = -1 * K_here * (diff_1 + diff_2) * grad_dot(uu, ww);
                        // considering overall sign change in weak form
                        double bwu = 1 * K_here * (diff_1 + diff_2) * grad_dot(uu, ww);

                        // double bwu = grad_dot(uu, ww);
                        val += bwu * w * J;
                    }
                }

                if (val != 0) {
                    int i = linear_index(a, x, y) + 1;
                    int j = linear_index(b, x, y) + 1;
                    vals_buf.push_back({i, j, val});
                    // executor.synchronized([&]() { problem.add(i, j, val); });
                }
            }

            // 1's for Dirichlet BC
            for_boundary_dofs(x, y, [&](index_type dof) {
                if (is_fixed(dof, x, y)) {
                    int i = linear_index(dof, x, y) + 1;
                    vals_buf.push_back({i, i, 1});
                    // executor.synchronized([&]() { problem.add(i, i, 1); });
                }
            });

            executor.synchronized([&]() {
                for (auto& [i, j, val] : vals_buf) {
                    problem.add(i, j, val);
                }
            });
        });
    }

    bool is_fixed(index_type dof, const dimension& /*x*/, const dimension& /*y*/) const {
        // return false; //dof[0] == 0 && dof[1] == 0; //|| dof[1] == y.dofs() - 1;
        return dof[0] == 0 || dof[0] == x.dofs() - 1 || dof[1] == 0 || dof[1] == y.dofs() - 1;
    }

    // bool is_fixed(index_type dof, const dimension& /*x*/, const dimension& /*y*/) const {
    //     return dof[0] == 0 || dof[0] == x.dofs() - 1 || dof[1] == 0 || dof[1] == y.dofs() - 1;
    // }
};

}  // namespace ads::problems

#endif  // CO2_SEQUESTRATION_CO2_SEQUESTRATION_2D_HPP
