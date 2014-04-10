/** Copyright (C) 2012 Benjamin Charrow
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#ifndef GMM_H
#define GMM_H

#include <cstdio>
#include <cmath>              
#include <algorithm>
#include <vector>
#include <list>

#include <ros/ros.h>
#include <ros/assert.h>

#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/Cholesky>
#include <Eigen/StdVector>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_integration.h>

#include <timeit/timeit.hpp>

#include <ros_gsl/random.hpp>

namespace rf {
  /* Approximate, but fast, exponentiation.
     http://martin.ankerl.com/2007/10/04/optimized-pow-approximation-for-java-and-c-c
  */
  double fast_exp(double val);
  
  typedef Eigen::Matrix<double, 1, 1> Matrix1d;
  typedef Eigen::Matrix<double, 1, 1> Vector1d; 

  template <int D>
  class GaussMM {
  public:
    // Types
    typedef Eigen::Matrix<double, D, 1> VectorD;
    typedef Eigen::Matrix<double, D, D> MatrixD;
    typedef Eigen::DiagonalWrapper<VectorD> CovWrap;
    typedef std::vector<VectorD> VectorList;
    typedef GaussMM<D>* Ptr;
    typedef const GaussMM<D>* ConstPtr;

    // Accessors
    const VectorList& means() const { return means_; }
    const VectorD& means(int i) const { return means_[i]; }
    const VectorList& covs() const { return covs_; }
    const VectorD& covs(int i) const { return covs_[i]; }    
    const VectorList& invCovs() const { return invcovs_; }
    const VectorD& invCovs(int i) const { return invcovs_[i]; }    
    const std::vector<double>& weights() const { return weights_; }
    double weights(int i) const { return weights_[i]; }
    
    size_t numComponents() const { return means_.size(); }
    static const int Dim = D;
    // This method is overriden when D=Eigen::Dynamic
    int dim() const {
      return Dim;
    }

    // Constructors
    GaussMM();
    GaussMM(const GaussMM &other);
    ~GaussMM();

    // Likelihoods
    double compLikelihood(int i) const {
      return normalizers_[i];
    }
    
    double compLikelihood(const VectorD &vec, int i) const {
      VectorD diff = vec - means_[i];
      double quad = diff.cwiseProduct(invcovs_[i]).cwiseProduct(diff).sum();
      return normalizers_[i] * exp(-0.5 * quad);
    }
    
    double likelihood(const VectorD &vec) const {
      double likelihood = 0.0;
      for (size_t i = 0; i < numComponents(); ++i) {
        VectorD diff = vec - means_[i];
        double quad = diff.cwiseProduct(invcovs_[i]).cwiseProduct(diff).sum();
        likelihood += normalizers_[i] * exp(-0.5 * quad);
      }
      return likelihood;
    };

    // Factory methods
    template<int NewD>
    static typename GaussMM<NewD>::Ptr Combine(const std::vector<GaussMM<D>::Ptr> &mixtures);
    static Ptr Sum(const std::vector<GaussMM<D>::Ptr> &mixtures,
                   const std::vector<double> &weights);

    // Mutators
    void addComponent(const VectorD &mean, const VectorD &cov, double weight) {
      means_.push_back(mean);
      covs_.push_back(cov);
      invcovs_.push_back(cov.cwiseInverse());
      weights_.push_back(weight);
      double det = cov.prod();
      normalizers_.push_back(weight / (pow(2 * M_PI, mean.size() / 2.0) * sqrt(det)));
    }

    void splitComponent(int i);
    Eigen::VectorXd sample(int num) const;
  private:
    static void findCombos(const std::vector<GaussMM::Ptr> &mixtures, size_t index,
                           std::vector<int> *vec, std::list<std::vector<int> > *combos);

    VectorList means_, covs_, invcovs_;
    std::vector<double> weights_, normalizers_;
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF_VECTORIZABLE_FIXED_SIZE(double, D);
  };

  typedef GaussMM<1> GaussMM1;
  typedef GaussMM<2> GaussMM2;

  //======================= Template specializations ========================//
  template<> int GaussMM<Eigen::Dynamic>::dim() const;
  template<> double GaussMM<1>::compLikelihood(const Vector1d &vec, int i) const;
  template<> double GaussMM<2>::compLikelihood(const Eigen::Vector2d &vec, int i) const;
  template <> Eigen::VectorXd GaussMM<1>::sample(int num) const;
  
  //========================== Method definitions ===========================//
  template <int D>
  std::ostream& operator<<(std::ostream &os, const GaussMM<D> &gmm) {
    for (size_t i = 0; i < gmm.numComponents(); ++i) {
      os << "Component " << i << std::endl;
      os << "  Weight: " << gmm.weights()[i] << std::endl;
      os << "  Mean: " << gmm.means(i).transpose() << std::endl;
      os << "  Cov: " << gmm.covs(i).transpose() << std::endl;
    }
    return os;
  }

  template <int D>
  GaussMM<D>::GaussMM() {
  
  }

  template <int D>
  GaussMM<D>::GaussMM(const GaussMM<D> &other) {
    means_ = other.means_;
    covs_ = other.covs_;
    invcovs_ = other.invcovs_;
    weights_ = other.weights_;
    normalizers_ = other.normalizers_;
  }

  template <int D>
  GaussMM<D>::~GaussMM() {
    
  }

  template <int D>
  template <int NewD>
  typename GaussMM<NewD>::Ptr GaussMM<D>::Combine(const std::vector<GaussMM<D>::Ptr> &mixtures) {

    std::list<std::vector<int> > combos;
    findCombos(mixtures, 0, NULL, &combos);

    // Debugging
    // fprintf(stderr, "Total lists: %i\n", combos.size());
    // for (std::list<std::vector<int> >::iterator it = combos.begin();
    //      it != combos.end();
    //      ++it) {
    //   std::vector<int> vec = (*it);
    //   for (std::vector<int>::const_iterator it = vec.begin(); it != vec.end(); ++it) {
    //     fprintf(stderr, "% 2i ", *it);
    //   }
    //   fprintf(stderr, "\n");
    // }
  
    // Build the new mixture model!
    typename GaussMM<NewD>::Ptr gmm = new GaussMM<NewD>();
    for (std::list<std::vector<int > >::iterator it = combos.begin();
         it != combos.end();
         ++it) {
      Eigen::Matrix<double, NewD, 1> mean, cov;

      if (NewD != Eigen::Dynamic) {
        mean = Eigen::Matrix<double, NewD, 1>::Zero(NewD);
        cov = Eigen::Matrix<double, NewD, 1>::Zero(NewD);
      } else {
        mean = Eigen::Matrix<double, NewD, 1>::Zero(D * mixtures.size());
        cov = Eigen::Matrix<double, NewD, 1>::Zero(D * mixtures.size());
      }
      
      double weight = 1.0;
      std::vector<int> indices = *it;
    
      // Construct new element out of components from previous mixtures
      int mat_index = 0;
      for (size_t i = 0; i < indices.size(); ++i) {
        typename GaussMM<D>::VectorD comp_mean = mixtures[i]->means(indices[i]);
        typename GaussMM<D>::VectorD comp_cov = mixtures[i]->covs(indices[i]);
        double comp_weight = mixtures[i]->weights(indices[i]);
        mean.segment(mat_index, D) = comp_mean;
        cov.segment(mat_index, D) = comp_cov;
        weight *= comp_weight;
      
        mat_index += D;
      }
    
      gmm->addComponent(mean, cov, weight);
    }
    return gmm;
  }

  template <int D>
  void GaussMM<D>::findCombos(const std::vector<Ptr> &mixtures, size_t index,
                              std::vector<int> *vec, std::list<std::vector<int> > *combos) {
    if (index == 0) {
      vec = new std::vector<int>();
    }
  
    if (index >= mixtures.size()) {
      combos->push_back(*vec);
      return;
    }

    ConstPtr mixture = mixtures.at(index);
    for (size_t i = 0; i < mixture->numComponents(); ++i) {
      vec->push_back(i);
      findCombos(mixtures, index + 1, vec, combos);
      vec->pop_back();
    }

    if (index == 0) {
      delete vec;
    }
  }

  template <int D>
  typename GaussMM<D>::Ptr GaussMM<D>::Sum(const std::vector<Ptr> &mixtures,
                                           const std::vector<double> &weights) {
    if (mixtures.size() != weights.size()) {
      ROS_ERROR("GaussMM::Sum() Weights and mixtures are different lengths");
    }
    Ptr gmm = new GaussMM();
  
    typename std::vector<Ptr>::const_iterator mix_it = mixtures.begin();
    std::vector<double>::const_iterator weight_it = weights.begin();
    for (; mix_it != mixtures.end() && weight_it != weights.end();
         ++weight_it, ++mix_it) {
      double weight = *weight_it;
      typename GaussMM<D>::Ptr mix = *mix_it;
      for (size_t i = 0; i < mix->numComponents(); ++i) {
        gmm->addComponent(mix->means(i), mix->covs(i), mix->weights(i) * weight);
      }
    }
    return gmm;
  }

  template <int D>
  void GaussMM<D>::splitComponent(int comp_ind) {
    double weights[4] = {0.12738084098, 0.37261915901, 0.37261915901, 0.12738084098};
    double means[4] = {-1.4131205233, -0.44973059608, 0.44973059608, 1.4131205233};
    double sd = 0.51751260421;
    
    // Find dimension with maximum variance
    VectorD eigs = covs_[comp_ind];
    double max_eig = -1;
    int max_ind = -1;
    for (int i = 0; i < dim(); ++i) {
      double eig = eigs(i);
      if (eig > max_eig) {
        max_eig = eig;
        max_ind = i;
      }
    }
    eigs(max_ind) *= sd * sd;

    const double orig_weight = weights_[comp_ind];
    const VectorD orig_mean = means_[comp_ind];
    VectorD e_d = VectorD::Zero(dim());
    e_d(max_ind) = 1;


    double det = sqrt(eigs.prod());
    for (int i = 0; i < 4; ++i) {
      VectorD mean = orig_mean + (means[i] * sqrt(max_eig) * e_d);
      // Put first new component where original component was, rest get added
      // normally.  At end of loop, free original component
      if (i == 0) {
        means_[comp_ind] = mean;
        covs_[comp_ind] = eigs;
        invcovs_[comp_ind] = eigs.cwiseInverse();
        weights_[comp_ind] = weights[i] * orig_weight;
        normalizers_[comp_ind] = weights_[comp_ind] / (pow(2 * M_PI, mean.size() / 2.0) * det);
      } else {
        addComponent(mean, eigs, weights[i] * orig_weight);
      }
    }
  }
  
  //================================= Entropy =================================//

  template <int D>
  class EntropyApproximator {
  public:
    EntropyApproximator() { ;}
    EntropyApproximator(const EntropyApproximator<D> &ea) { ; }    
    virtual ~EntropyApproximator() { ; }

    virtual double evaluate(const GaussMM<D> &gmm) = 0;

    // static double gaussEntropy(Eigen::Matrix<double, D, D> cov);
    static double gaussEntropy(Eigen::Matrix<double, D, 1> cov);
  };

  template <int D>
  class NumericalApproximator : public EntropyApproximator<D> {
  public:
    NumericalApproximator(int iters) : iters_(iters) {
      ;
    }
    NumericalApproximator(const NumericalApproximator<D> &na) :
      EntropyApproximator<D>(na), iters_(na.iters_) {
      ;
    }
    
    double evaluate(const GaussMM<D> &gmm);
  private:
    int iters_;
  };

  
  template <int D>
  class TaylorApproximator : public EntropyApproximator<D> {
  public:
    TaylorApproximator(int order, int num_splits) :
      order_(order), splits_(num_splits) {
      if (order_ < 0 || order_ > 4) {
        ROS_ERROR("Order must be in range [0, 4]");
        ROS_BREAK();
      }
      if (splits_ < 0) {
        ROS_ERROR("Splits must be positive integer");
        ROS_BREAK();
      }
    }

    TaylorApproximator(const TaylorApproximator<D> &ta) :
      EntropyApproximator<D>(ta), order_(ta.order_), splits_(ta.splits_) {
      ;
    }
    
    double evaluate(const GaussMM<D> &gmm);
  private:
    typedef std::vector<typename GaussMM<D>::VectorD,
                        Eigen::aligned_allocator<typename GaussMM<D>::VectorD> > EigStdVector;

    int order_, splits_;
    
    double taylorInternal(typename GaussMM<D>::ConstPtr gmm);
    EigStdVector gradientAtMeans(typename GaussMM<D>::ConstPtr gmm);
  };

  template <int D>
  double EntropyApproximator<D>::gaussEntropy(Eigen::Matrix<double, D, 1> cov) {
    return log(pow(2 * M_PI * exp(1), cov.size() / 2.0) * sqrt(cov.prod()));
  }
  
  template <int D>
  double eval_gmm_ent(double z[], size_t dim, void *g) {
    typename GaussMM<D>::ConstPtr gmm = static_cast<typename GaussMM<D>::ConstPtr>(g);
    typename GaussMM<D>::VectorD vec(dim);
    for (size_t i = 0; i < dim; ++i) {
      vec(i) = z[i];
    }
    double likelihood = gmm->likelihood(vec);
    if (likelihood < 1e-50) {
      return 0;
    } else {
      return -likelihood * log(likelihood);
    }
  }
  
    
  template <int D>
  double NumericalApproximator<D>::evaluate(const GaussMM<D> &gmm) {
    // timeit::toc("Test!");
    double result;
    gsl_monte_function G;
    G.f= &eval_gmm_ent<D>;
    G.dim = gmm.dim();
    G.params = static_cast<void*>(const_cast<typename GaussMM<D>::Ptr>(&gmm));

    double *lower, *upper;
    lower = new double[gmm.dim()];
    upper = new double[gmm.dim()];
    for (int i = 0; i < gmm.dim(); ++i) {
      double mean_max = -std::numeric_limits<double>::infinity();
      double mean_min = std::numeric_limits<double>::infinity();
      double cov_max = -1;

      for (size_t comp = 0; comp < gmm.numComponents(); ++comp) {
        mean_max = std::max(mean_max, gmm.means(comp)(i));
        mean_min = std::min(mean_min, gmm.means(comp)(i));
        cov_max = std::max(cov_max, gmm.covs(comp)(i));
      }

      lower[i] = mean_min - 10 * sqrt(cov_max);
      upper[i] = mean_max + 10 * sqrt(cov_max);
    }
    
    gsl_monte_vegas_state *vegas;
    vegas = gsl_monte_vegas_alloc(gmm.dim());
    gsl_monte_vegas_init(vegas);

    double error;
    gsl_monte_vegas_integrate(&G, lower, upper, gmm.dim(), iters_,
                              gsl::rng(), vegas, &result, &error);
    gsl_monte_vegas_free(vegas);    
    // timeit::toc("NA:post_int!");
    return result;
  }

  template <> double NumericalApproximator<1>::evaluate(const GaussMM<1> &gmm);
  
  template <int D>
  double TaylorApproximator<D>::evaluate(const GaussMM<D> &gmm) {
    // timeit::toc("TA:start");
    if (splits_ == 0) {
      return taylorInternal(&gmm);
    } else {
      // Split stored GMM splits times
      typename GaussMM<D>::Ptr g = new GaussMM<D>(gmm);
      for (int i = 0; i < splits_; ++i) {
        // Find component with biggest variance
        int max_var_ind = -1;
        double max_var = -1.0;
        const typename GaussMM<D>::VectorList & covs = g->covs();
        
        for (size_t j = 0; j < covs.size(); ++j) {
          double var = covs[j].maxCoeff();
          if (max_var < var) {
            max_var_ind = j;
            max_var = var;
          }
        }
        g->splitComponent(max_var_ind);
      }
      double result = taylorInternal(g);
      delete g;
      return result;
    }
  }

  template <int D>
  double TaylorApproximator<D>::taylorInternal(typename GaussMM<D>::ConstPtr g) {
    // timeit::tic();
    // timeit::toc("TA:taylor!");
    double entropy = 0.0;
    // Evaluate logarithim of likelihood at mean of each component
    std::vector<double> gmm_likelihoods;
    gmm_likelihoods.resize(g->numComponents());

    const std::vector<double>& weights = g->weights();
    const typename GaussMM<D>::VectorList& means = g->means();
    const typename GaussMM<D>::VectorList& covs = g->covs();
    const typename GaussMM<D>::VectorList& invcovs = g->invCovs();
    // timeit::toc("TA:taylor1!");  
    for (size_t ind = 0; ind < weights.size(); ++ind) {
      gmm_likelihoods[ind] = g->likelihood(means[ind]);
      // Compute entropy for 0 order expansion
      entropy += -weights[ind] * log(gmm_likelihoods[ind]);
    }
    // timeit::toc("TA:taylor3!");
    // 2nd order expansion
    double entropy2 = 0.0;
    if (order_ >= 2) {
      EigStdVector grads = gradientAtMeans(g);
      // timeit::toc("TA:grads");
      for (size_t k = 0; k < g->numComponents(); ++k) {
        typename GaussMM<D>::VectorD mean = means[k];
        double mean_likelihood = gmm_likelihoods[k];
      
        typename GaussMM<D>::VectorD total = GaussMM<D>::VectorD::Zero(g->dim());

        for (size_t j = 0; j < g->numComponents(); ++j) {
          typename GaussMM<D>::VectorD scale = invcovs[j];
          typename GaussMM<D>::VectorD mu_diff = mean - means[j];
          typename GaussMM<D>::VectorD mu_mu_outer = mu_diff.cwiseProduct(mu_diff);
          typename GaussMM<D>::VectorD mu_grad_outer = mu_diff.cwiseProduct(grads[k]);  

          typename GaussMM<D>::VectorD term = (GaussMM<D>::VectorD::Ones(g->dim()) -
                                               mu_mu_outer.cwiseProduct(scale) -
                                               mu_grad_outer / mean_likelihood);
          total -= scale.cwiseProduct(term) * g->compLikelihood(mean, j);
        }
      
        double sum = total.cwiseProduct(covs[k]).sum() / mean_likelihood;
        entropy2 -= 0.5 * weights[k] * sum;
      }
    }
    // timeit::toc("TA:done!");  
    return entropy + entropy2;
  }

  template <int D>
  typename TaylorApproximator<D>::EigStdVector
  TaylorApproximator<D>::gradientAtMeans(typename GaussMM<D>::ConstPtr g) {
    EigStdVector grads;
    grads.resize(g->numComponents());
    for (size_t i = 0; i < g->numComponents(); ++i) {
      grads[i].resize(D);
      grads[i].setZero();
    }

    const typename GaussMM<D>::VectorList& means = g->means();
    const typename GaussMM<D>::VectorList& invcovs = g->invCovs();
    
    for (size_t j = 0; j < g->numComponents(); ++j) {
      typename GaussMM<D>::VectorD mean = means[j];
      typename GaussMM<D>::VectorD invcov = invcovs[j];
      for (size_t i = 0; i < g->numComponents(); ++i) {
        typename GaussMM<D>::VectorD data = means[i];
        grads[i] -= (data - mean).cwiseProduct(invcov) * g->compLikelihood(data, j);
      }
    }
    return grads;
  }
}
#endif
