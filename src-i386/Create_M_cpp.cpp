#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


// //' @export
// [[Rcpp::export]]
NumericVector CreateM_cpp(IntegerVector J_s, IntegerVector Number_of_Est_s, IntegerVector K_s, NumericVector g_s, NumericVector epsilon_s, NumericVector m_s, NumericVector s_squared_s, NumericVector v_s, int n) {
  int i;
  int j;
  int k;
  int l;
  int i_1;
  int i_2;
  int j_1;
  int j_2;
  int k_1;
  int k_2;
  int mu;
  int nu;
  int Counter = 0;

  IntegerVector Number_of_Est (n+1);
  IntegerVector Number_of_Est_cum (n+1);

  int Max_Number_of_Est = 0;
  Number_of_Est_cum[0] = 0;
  for (i = 1; i <= n; i++) {
    Number_of_Est[i] = Number_of_Est_s[i-1];
    Number_of_Est_cum[i] =  Number_of_Est_cum[i-1] + Number_of_Est[i];
    if (Number_of_Est[i] > Max_Number_of_Est) {
      Max_Number_of_Est = Number_of_Est[i];
    }
  }

  IntegerMatrix J_all (Number_of_Est_cum[n]+1, 4);
  for (i=1; i<=Number_of_Est_cum[n]; i++) {
    J_all(i,1) = J_s[i-1];
    J_all(i,2) = J_s[i-1+Number_of_Est_cum[n]];
    J_all(i,3) = J_s[i-1+2*Number_of_Est_cum[n]];
  }

  int Count_K = K_s.size();
  IntegerVector K (Count_K + 1);
  for (k=1; k <= Count_K; k++) {
    K[k] = K_s[k-1];
  }


  NumericMatrix g (n+1, Max_Number_of_Est + 1); // g(i,j) entspricht g[[i]][j]
  for (i = 1; i <= n; i++) {
    for (j = Number_of_Est_cum[i-1]; j < Number_of_Est_cum[i]; j++) {
      g(i, j+1-Number_of_Est_cum[i-1]) = g_s[j];
    }
  }

  NumericVector epsilon (n);
  for (i=1; i<=n-1; i++) {
    epsilon[i] = epsilon_s[i-1];
  }
  NumericVector m (n+1);
  for (k=1; k<=n; k++) {
    m[k] = m_s[k-1];
  }


  NumericVector s_squared (n+1);
  for (k=1; k<=n; k++) {
    s_squared[k] = s_squared_s[k-1];
  }
  NumericVector v (n+1);
  for (i=1; i<=n; i++) {
    v[i] = v_s[i-1];
  }


  // calculate gammas
  // set gam[i,l] <- 1 if $l\not\in K$

  NumericMatrix gam (n + 1, n + 1);
  for (i = 1; i <= n; i++) {
    for (k = 1; k <= n; k++) {
      gam(i,k) = 1;
    }
  }
  for (i = 1; i <= n; i++) {
    for (k = 1; k <= Count_K; k++) {
      gam(i,K[k]) = 1+s_squared[K[k]]/(v[i]*m[K[k]]*m[K[k]]);
    }
  }

  // calculate squared epsilons
  NumericVector epsilon_sqared (n);
  for (i=1; i<=n-1; i++) {
    epsilon_sqared[i] = epsilon[i] * epsilon[i];
  }

  // calculate \prod_{\nu=i}^{n-1} (1+\varepsilon_\nu^2)
  // set to 1 if i=n
  NumericVector prod1plusepssq (n+1);
  prod1plusepssq (n) = 1;
  for (i=n-1; i>=1; i--) {
    prod1plusepssq[i] = prod1plusepssq[i+1] * (1+epsilon_sqared[i]);
  }

  k=0;
  for (l=1; l<=n; l++) {
    k += (n-l+1)*(n-l+1);
  }
  NumericVector M_vector(k);

  for (l=1; l<=n; l++) {
    NumericMatrix CovMatrix(Number_of_Est_cum[n-l+1]+1, Number_of_Est_cum[n-l+1]+1);

    for (nu =1; nu <= Number_of_Est_cum[n-l+1]; nu++) {
      for (mu =nu; mu <=Number_of_Est_cum[n-l+1]; mu++) {
        i_1 = J_all(nu, 1);
        i_2 = J_all(mu, 1);
        j_1 = J_all(nu, 2);
        j_2 = J_all(mu, 2);
        k_1 = J_all(nu, 3);
        k_2 = J_all(mu, 3);
        if (i_1 == i_2)  {
          if (j_1 == 0 && j_2 == 0) {
            CovMatrix(nu,mu) = (m[l]*m[l] + s_squared[l]/v[i_1]) * prod1plusepssq[i_1] - m[l]*m[l];
          } else if (j_1 == 0 && k_2 == l) {
            CovMatrix(nu,mu) = m[l]*m[l]/gam(i_1,l) * (prod1plusepssq[std::max(i_1,j_2)]-1);
          } else if (j_1 == 0 && k_2 != l) {
            CovMatrix(nu,mu) = (m[l]*m[l] + s_squared[l]/v[i_1]) * prod1plusepssq[std::max(i_1,j_2)] - m[l]*m[l];
          } else if (j_2 == 0 && k_1 == l) {
            CovMatrix(nu,mu) = m[l]*m[l]/gam(i_1,l) * (prod1plusepssq[std::max(i_1,j_1)]-1);
          } else if (j_2 == 0 && k_1 != l) {
            CovMatrix(nu,mu) = (m[l]*m[l] + s_squared[l]/v[i_1]) * prod1plusepssq[std::max(i_1,j_1)] - m[l]*m[l];
          } else if (j_1 == j_2 && k_1 == k_2 && k_1 == l) {
            CovMatrix(nu,mu) = m[l]*m[l]/(gam(i_1,l)*gam(i_1,l)) * (gam(j_1,l)*prod1plusepssq[j_1]-1);
          } else if (j_1 == j_2 && k_1 == k_2 && k_1 != l) {
            CovMatrix(nu,mu) = (m[l]*m[l] + s_squared[l]/v[i_1]) * gam(i_1,k_1)*gam(j_1,k_1)*prod1plusepssq[j_1] - m[l]*m[l];
          } else if (j_1 != j_2 && k_1 == k_2 && k_1 != l) {
            CovMatrix(nu,mu) = (m[l]*m[l] + s_squared[l]/v[i_1]) * gam(i_1,k_1)*prod1plusepssq[std::max(j_1,j_2)] - m[l]*m[l];
          } else if (j_1 != j_2 && k_1 == k_2 && k_1 == l) {
            CovMatrix(nu,mu) = m[l]*m[l]/(gam(i_1,l)*gam(i_1,l)) * (prod1plusepssq[std::max(j_1,j_2)]-1);
          } else if (k_1 != k_2 && k_1 != l && k_2 != l) {
            CovMatrix(nu,mu) = (m[l]*m[l]+s_squared[l]/v[i_1]) * prod1plusepssq[std::max(j_1,j_2)] - m[l]*m[l];
          } else if (k_1 != k_2 && k_1 == l) {
            CovMatrix(nu,mu) = m[l]*m[l]/gam(i_1,l) * (prod1plusepssq[std::max(j_1,j_2)]-1);
          } else if (k_1 != k_2 && k_2 == l) {
            CovMatrix(nu,mu) = m[l]*m[l]/gam(i_1,l) * (prod1plusepssq[std::max(j_1,j_2)]-1);
          } else  {
            Rcout << "Fehler1";
          }
        } else {
          if (j_1 == 0 && j_2 == 0) {
            CovMatrix(nu,mu) = m[l]*m[l] * (prod1plusepssq[std::max(i_1,i_2)]-1);
          } else if (j_2 == 0 && j_1 == i_2 && k_1 == l) {
            CovMatrix(nu,mu) = m[l]*m[l]/gam(i_1,l) * (gam(i_2,l) * prod1plusepssq[i_2]-1);
          } else if (j_1 == 0 && j_2 == i_1 && k_2 == l) {
            CovMatrix(nu,mu) = m[l]*m[l]/gam(i_2,l) * (gam(i_1,l) * prod1plusepssq[i_1]-1);
          } else if (j_2 == 0 && j_1 != i_2 && k_1 == l) {
            CovMatrix(nu,mu) = m[l]*m[l]/gam(i_1,l) * (prod1plusepssq[std::max(j_1,i_2)]-1);
          } else if (j_1 == 0 && j_2 != i_1 && k_2 == l) {
            CovMatrix(nu,mu) = m[l]*m[l]/gam(i_2,l) * (prod1plusepssq[std::max(j_2,i_1)]-1);
          } else if (j_2 == 0 && k_1 != l) {
            CovMatrix(nu,mu) = m[l]*m[l] * (prod1plusepssq[std::max(j_1,i_2)]-1);
          } else if (j_1 == 0 && k_2 != l) {
            CovMatrix(nu,mu) = m[l]*m[l] * (prod1plusepssq[std::max(j_2,i_1)]-1);
          } else if (j_1 == j_2 && k_1 == k_2 && k_1 == l) {
            CovMatrix(nu,mu) = m[l]*m[l]/(gam(i_1,l)*gam(i_2,l)) * (gam(j_1,l) * prod1plusepssq[j_1]-1);
          } else if (j_1 == j_2 && k_1 == k_2 && k_1 != l) {
           CovMatrix(nu,mu) = m[l]*m[l] * (gam(j_1,k_1) * prod1plusepssq[j_1]-1);
          } else if (j_1 != j_2 && j_1 != i_2 && j_2 != i_1 && k_1 == k_2 && k_1 != l) {
            CovMatrix(nu,mu) = m[l]*m[l] * (prod1plusepssq[std::max(j_1,j_2)]-1);
          } else if (j_1 != j_2 && j_1 == i_2 && j_2 == i_1 && k_1 == k_2 && k_1 != l) {
            CovMatrix(nu,mu) = m[l]*m[l] * (prod1plusepssq[std::max(j_1,j_2)]/(gam(i_1,k_1)*gam(i_2,k_2))-1);
          } else if (j_1 != j_2 && j_1 != i_2 && j_2 == i_1 && k_1 == k_2 && k_1 != l) {
            CovMatrix(nu,mu) = m[l]*m[l] * (prod1plusepssq[std::max(j_1,j_2)]/gam(i_1,k_1)-1);
          } else if (j_1 != j_2 && j_1 == i_2 && j_2 != i_1 && k_1 == k_2 && k_1 != l) {
            CovMatrix(nu,mu) = m[l]*m[l] * (prod1plusepssq[std::max(j_1,j_2)]/gam(i_2,k_2)-1);
          } else if (j_1 != j_2 && k_1 == k_2 && k_1 == l) {
            CovMatrix(nu,mu) = m[l]*m[l]/(gam(i_1,l)*gam(i_2,l)) * (prod1plusepssq[std::max(j_1,j_2)]-1);
          } else if (j_1 != i_2 && k_1 == l && k_2 != l) {
            CovMatrix(nu,mu) = m[l]*m[l]/gam(i_1,l) * (prod1plusepssq[std::max(j_1,j_2)]-1);
          } else if (j_2 != i_1 && k_2 == l && k_1 != l) {
            CovMatrix(nu,mu) = m[l]*m[l]/gam(i_2,l) * (prod1plusepssq[std::max(j_1,j_2)]-1);
          } else if (j_1 == i_2 && k_1 == l && k_2 != l) {
            CovMatrix(nu,mu) = m[l]*m[l]/gam(i_1,l) * (gam(i_2,l)*prod1plusepssq[std::max(j_1,j_2)]-1);
          } else if (j_2 == i_1 && k_2 == l && k_1 != l) {
            CovMatrix(nu,mu) = m[l]*m[l]/gam(i_2,l) * (gam(i_1,l)*prod1plusepssq[std::max(j_1,j_2)]-1);
          } else if (k_1 != l && k_2 != l && k_1 != k_2) {
            CovMatrix(nu,mu) = m[l]*m[l] * (prod1plusepssq[std::max(j_1,j_2)]-1);
          } else {
            Rcout << "Fehler2";
          }
        }
      }
    }
    for (nu=2; nu<=Number_of_Est_cum[n-l+1]; nu++) {
      for (mu = 1; mu <nu; mu++) {
        CovMatrix(nu, mu) = CovMatrix(mu, nu);
      }
    }

    NumericMatrix M_l(n-l+2, n-l+2);
    for (i=1; i<=n-l+1; i++) {
      for (j=i; j<=n-l+1; j++) {
        M_l(i,j) = 0;
        for (nu=1; nu <= Number_of_Est[i]; nu++) {
          for (mu=1; mu <=Number_of_Est[j]; mu++) {
            M_l(i,j) = M_l(i,j) + g(i,nu) * CovMatrix(Number_of_Est_cum[i-1]+nu, Number_of_Est_cum[j-1]+mu) * g(j,mu);
          }
        }
        if (j != i) {
          M_l(j,i) = M_l(i,j);
        }
      }
    }

    for (i=1; i<=n-l+1; i++) {
      for (j=1; j<=n-l+1; j++) {
        M_vector[Counter] = M_l(i,j);
        Counter++;
      }
    }
  }

  return(M_vector);
}


