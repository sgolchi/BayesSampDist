// Generated by rstantools.  Do not edit by hand.

/*
    BayesSampDist is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    BayesSampDist is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with BayesSampDist.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.21.0
#include <stan/model/model_header.hpp>
namespace model_fitab_H0_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_fitab_H0");
    reader.add_event(18, 16, "end", "model_fitab_H0");
    return reader;
}
#include <stan_meta_header.hpp>
class model_fitab_H0
  : public stan::model::model_base_crtp<model_fitab_H0> {
private:
        int N;
        std::vector<double> a;
        std::vector<int> n;
public:
    model_fitab_H0(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_fitab_H0(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_fitab_H0_namespace::model_fitab_H0";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 2;
            context__.validate_dims("data initialization", "N", "int", context__.to_vec());
            N = int(0);
            vals_i__ = context__.vals_i("N");
            pos__ = 0;
            N = vals_i__[pos__++];
            check_greater_or_equal(function__, "N", N, 1);
            current_statement_begin__ = 3;
            validate_non_negative_index("a", "N", N);
            context__.validate_dims("data initialization", "a", "double", context__.to_vec(N));
            a = std::vector<double>(N, double(0));
            vals_r__ = context__.vals_r("a");
            pos__ = 0;
            size_t a_k_0_max__ = N;
            for (size_t k_0__ = 0; k_0__ < a_k_0_max__; ++k_0__) {
                a[k_0__] = vals_r__[pos__++];
            }
            size_t a_i_0_max__ = N;
            for (size_t i_0__ = 0; i_0__ < a_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "a[i_0__]", a[i_0__], 0);
                check_less_or_equal(function__, "a[i_0__]", a[i_0__], 1);
            }
            current_statement_begin__ = 5;
            validate_non_negative_index("n", "N", N);
            context__.validate_dims("data initialization", "n", "int", context__.to_vec(N));
            n = std::vector<int>(N, int(0));
            vals_i__ = context__.vals_i("n");
            pos__ = 0;
            size_t n_k_0_max__ = N;
            for (size_t k_0__ = 0; k_0__ < n_k_0_max__; ++k_0__) {
                n[k_0__] = vals_i__[pos__++];
            }
            size_t n_i_0_max__ = N;
            for (size_t i_0__ = 0; i_0__ < n_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "n[i_0__]", n[i_0__], 0);
            }
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 8;
            num_params_r__ += 1;
            current_statement_begin__ = 9;
            num_params_r__ += 1;
            current_statement_begin__ = 10;
            num_params_r__ += 1;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_fitab_H0() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 8;
        if (!(context__.contains_r("alpha1")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable alpha1 missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("alpha1");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "alpha1", "double", context__.to_vec());
        double alpha1(0);
        alpha1 = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, alpha1);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable alpha1: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 9;
        if (!(context__.contains_r("alpha2")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable alpha2 missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("alpha2");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "alpha2", "double", context__.to_vec());
        double alpha2(0);
        alpha2 = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, alpha2);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable alpha2: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 10;
        if (!(context__.contains_r("sigma_a")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable sigma_a missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("sigma_a");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "sigma_a", "double", context__.to_vec());
        double sigma_a(0);
        sigma_a = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, sigma_a);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable sigma_a: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 8;
            local_scalar_t__ alpha1;
            (void) alpha1;  // dummy to suppress unused var warning
            if (jacobian__)
                alpha1 = in__.scalar_lb_constrain(0, lp__);
            else
                alpha1 = in__.scalar_lb_constrain(0);
            current_statement_begin__ = 9;
            local_scalar_t__ alpha2;
            (void) alpha2;  // dummy to suppress unused var warning
            if (jacobian__)
                alpha2 = in__.scalar_lb_constrain(0, lp__);
            else
                alpha2 = in__.scalar_lb_constrain(0);
            current_statement_begin__ = 10;
            local_scalar_t__ sigma_a;
            (void) sigma_a;  // dummy to suppress unused var warning
            if (jacobian__)
                sigma_a = in__.scalar_lb_constrain(0, lp__);
            else
                sigma_a = in__.scalar_lb_constrain(0);
            // model body
            current_statement_begin__ = 13;
            for (int i = 1; i <= N; ++i) {
                current_statement_begin__ = 14;
                lp_accum__.add(lognormal_log<propto__>(get_base1(a, i, "a", 1), ((-(alpha1) / get_base1(n, i, "n", 1)) + (alpha2 / pow(get_base1(n, i, "n", 1), 2))), sigma_a));
            }
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("alpha1");
        names__.push_back("alpha2");
        names__.push_back("sigma_a");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_fitab_H0_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        double alpha1 = in__.scalar_lb_constrain(0);
        vars__.push_back(alpha1);
        double alpha2 = in__.scalar_lb_constrain(0);
        vars__.push_back(alpha2);
        double sigma_a = in__.scalar_lb_constrain(0);
        vars__.push_back(sigma_a);
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            if (!include_gqs__ && !include_tparams__) return;
            if (!include_gqs__) return;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    std::string model_name() const {
        return "model_fitab_H0";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "alpha1";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "alpha2";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "sigma_a";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "alpha1";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "alpha2";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "sigma_a";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
    }
}; // model
}  // namespace
typedef model_fitab_H0_namespace::model_fitab_H0 stan_model;
#ifndef USING_R
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
#endif
#endif
