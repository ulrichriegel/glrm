#include <Rcpp.h>
#include <math.h>
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






double Calculate_E_Sr_cpp(IntegerVector & Index_S, IntegerVector & Index_r, NumericVector & prod1plusepssq, NumericMatrix & gam, IntegerVector & K, NumericVector &  m, NumericVector &  s_squared, NumericVector &  v, int n) {
    IntegerMatrix chi (n+1, n+1);
    IntegerMatrix xi (n+1, n+1);
    IntegerMatrix rho (n+1, n+1);

    Rcpp::checkUserInterrupt();
    int Count_r;
    int Count_S;
    int i;
    int j;
    int k;
    int l;

    double dblTemp;

    for (i=1; i<=n; i++) {
      for (j=1; j<=n; j++) {
        chi(i,j) = 0;
        xi(i,j) = 0;
       }
    }

    if (Index_r.size() == 4) {
      Count_r = 1;
    } else if (Index_r.size() == 7) {
      Count_r = 2;
    } else if (Index_r.size() == 10) {
      Count_r = 3;
    } else {
      Count_r = 4;
    }
    if (Index_S.size() == 3) {
      Count_S = 1;
    } else {
      Count_S = 2;
    }

    //Rcout << "Count_r: " << Count_r << "; Count_S: " << Count_S << "\n";

    for (int nu = 1; nu <= Count_r; nu++) {
      i = Index_r[3 * (nu - 1) + 1];
      j = Index_r[3 * (nu - 1) + 2];
      k = Index_r[3 * (nu - 1) + 3];
      //Rcout << "Index_r: i: " << i << "; j: " << j << "; k: " << k << "\n";
      if (j != 0) {
        chi(j,k) = chi(j,k) + 1;
        chi(i,k) = chi(i,k) - 1;
        xi(i,k) = xi(i,k) - 1;
      }
    }

    for (int nu = 1; nu <= Count_S; nu++) {
      i = Index_S[2 * (nu - 1) + 1];
      k = Index_S[2 * (nu - 1) + 2];
      //Rcout << "Index_S: i: " << i << "; k: " << k << "\n";
      for (j = 1; j <= K.size()-1; j++) {
        if (k == K[j]) {
          chi(i,k) = chi(i,k) + 1;
        }
      }
    }

    for (i=1; i<=n; i++) {
      for (j=1; j<=n; j++) {
        rho(i,j) = 0.5 * chi(i,j) * (chi(i,j) - 1);
      }
    }


    //Rcout << chi;
    //Rcout << "chi\n";
    //for (i = 1; i <= n; i++) {
    //  for (k = 1; k <= n; k++) {
    //    Rcout << chi(i,k) << " ";
    //  }
    //  Rcout << "\n";
    //}
    //Rcout << gam;
    //Rcout << "gam\n";
    //for (i = 1; i <= n; i++) {
    //  for (k = 1; k <= n; k++) {
    //    Rcout << gam(i,k) << " ";
    //  }
    //  Rcout << "\n";
    //}


    dblTemp = 1;
    for (i=1; i<=n; i++) {
      for (j=1; j<=n-1; j++) {
        if (rho(i,j)+xi(i,j) != 0) {
          dblTemp = dblTemp * std::pow(gam(i,j) , rho(i,j) + xi(i,j));
        }
      }
    }
//      Temp <- prod(gam^(rho+xi))

    IntegerVector pi (Count_r+1);
    for (i=1; i<=Count_r; i++) {
      if (Index_r[3*(i-1)+2]==0) {
        pi[i] = Index_r[3*(i-1)+1];
      } else {
        pi[i] = Index_r[3*(i-1)+2];
      }
    }
//      pi <- ifelse(Index_r[,2]==0,Index_r[,1],Index_r[,2])

    double dblBeta = 1;
    if (Count_r>1) {
      IntegerVector i_vec = clone(pi);
      i_vec[0] = -1;
      std::sort(i_vec.begin(), i_vec.end());
      for (int kap=2; kap<= Count_r; kap++) {
        dblBeta = dblBeta * std::pow(prod1plusepssq[i_vec[kap]], kap-1);
      }
    }
    //Rcout << dblBeta;
//      # Calculate beta
//      beta <- 1
//      if (Count_r > 1) {
//          std::sort(pi.begin(), pi.end());
//        i_vec <- sort(pi)
//        for (kappa in 2:Count_r) {
//          beta <- beta * prod1plusepssq[i_vec[kappa]]^(kappa-1)
//        }
//      }

    IntegerVector tau (Count_r+1);
    for (i=1; i<=Count_r; i++) {
      tau[i] = Index_r[3*(i-1)+1];
    }

    //Rcout << "Beta " << dblBeta;

//      tau <- Index_r[,1]
    bool blnInK = false;
    double E_Sr = dblTemp * dblBeta;
    for (i=1; i<=Count_r; i++) {
      E_Sr = E_Sr / v[tau[i]];
    }
    for (i=1; i<=Count_S; i++) {
      E_Sr = E_Sr * v[Index_S[2*(i-1)+1]];
    }

    if (Count_S > 1) {
      if (Index_S[1] == Index_S[3] && Index_S[2] == Index_S[4]) {
        for (i=1; i <= K.size()-1; i++) {
          if (Index_S[2] == K[i]) {
            blnInK = true;
          }
        }
        if (blnInK == false) {
          i = Index_S[1];
          l = Index_S[2];
          E_Sr = E_Sr * (m[l] * m[l] + s_squared[l] / v[i]);
        } else {
          for (i=1; i <= Count_S; i++) {
            E_Sr = E_Sr * m[Index_S[2*i]];
          }
          //
        }
      } else {
        for (i=1; i <= Count_S; i++) {
          E_Sr = E_Sr * m[Index_S[2*i]];
        }

      }
    } else {
      E_Sr = E_Sr * m[Index_S[2]];

    }
//      E_Sr <- beta / prod(v[tau]) * prod(v[Index_S[,1]]) * Temp
//        if (Count_S > 1) {
//          if ((Index_S[1,1] == Index_S[2,1]) && (Index_S[1,2] == Index_S[2,2]) && (Index_S[1,2] %in% (1:n)[-K])) {
//            i <- Index_S[1,1]
//            l <- Index_S[1,2]
//            E_Sr <- E_Sr * (m[l]^2 + s_squared[l] / v[i])
//          } else {
//            E_Sr <- E_Sr * prod(m[Index_S[,2]])
//          }
//        } else {
//          E_Sr <- E_Sr * prod(m[Index_S[,2]])
//        }
//        return(E_Sr)

  return (E_Sr);
}


double Cov_S_by_vhat_S_by_vhat_cpp(IntegerMatrix & J_all, NumericMatrix & g, IntegerVector & Index_S1, IntegerVector & Index_prodv1, int Dim1, IntegerVector & Index_S2, IntegerVector & Index_prodv2, int Dim2, NumericVector & prod1plusepssq, NumericMatrix & gam, IntegerVector & K, NumericVector & m, NumericVector & s_squared, NumericVector & v, int n, IntegerVector & Number_of_Est, IntegerVector & Number_of_Est_cum) {
  double Covariance = 0;
  double E_1;
  double E_2;
  double E_12;
  IntegerVector Index_S12 (5);
  for (int i=1; i<=2;i++) {
   Index_S12[i] = Index_S1[i];
   Index_S12[i+2] = Index_S2[i];
  }
  //int counter = 0;

  if (Dim1 == 1 && Dim2 == 1) {
    IntegerVector Index_r1 (4);
    IntegerVector Index_r2 (4);
    IntegerVector Index_r12 (7);
    for (int j1 = 1; j1 <= Number_of_Est[Index_prodv1[1]]; j1++) {
      Index_r1[1] = J_all (Number_of_Est_cum[Index_prodv1[1]-1]+j1,1);
      Index_r1[2] = J_all (Number_of_Est_cum[Index_prodv1[1]-1]+j1,2);
      Index_r1[3] = J_all (Number_of_Est_cum[Index_prodv1[1]-1]+j1,3);
      for (int k1 = 1; k1 <= Number_of_Est[Index_prodv2[1]]; k1++) {
        Index_r2[1] = J_all (Number_of_Est_cum[Index_prodv2[1]-1]+k1,1);
        Index_r2[2] = J_all (Number_of_Est_cum[Index_prodv2[1]-1]+k1,2);
        Index_r2[3] = J_all (Number_of_Est_cum[Index_prodv2[1]-1]+k1,3);
        for (int i=1; i<=3;i++) {
          Index_r12[i] = Index_r1[i];
          Index_r12[i+3] = Index_r2[i];
        }
        E_1 = Calculate_E_Sr_cpp(Index_S1, Index_r1, prod1plusepssq, gam, K, m, s_squared, v, n);
        E_2 = Calculate_E_Sr_cpp(Index_S2, Index_r2, prod1plusepssq, gam, K, m, s_squared, v, n);
        E_12 = Calculate_E_Sr_cpp(Index_S12, Index_r12, prod1plusepssq, gam, K, m, s_squared, v, n);
        //Rcout << "Index_S_12 " << Index_S12 << " Index_r12 " << Index_r12 << "\n";
        //counter++;
        //if (counter < 0) {
          //Rcout << "E1 " << E_1 << " E2 " << E_2 << " E12 " << E_12 << "\n";
          //Rcout << "g1 " << g(Index_prodv1[1], j1) << " g2 " << g(Index_prodv2[1], k1) << "\n";
        //}
        Covariance = Covariance + g(Index_prodv1[1], j1) * g(Index_prodv2[1], k1) * (E_12 - E_1 * E_2);
        //Rcout << "Cov " << counter << " " << Covariance << "\n";
      }
    }
  } else if (Dim1 == 2 && Dim2 == 2) {
    IntegerVector Index_r1 (7);
    IntegerVector Index_r2 (7);
    IntegerVector Index_r12 (13);
    for (int j1 = 1; j1 <= Number_of_Est[Index_prodv1[1]]; j1++) {
      for (int j2 = 1; j2 <= Number_of_Est[Index_prodv1[2]]; j2++) {
        for (int i=1; i<=3;i++) {
          Index_r1[i] = J_all (Number_of_Est_cum[Index_prodv1[1]-1]+j1,i);
          Index_r1[i+3] = J_all (Number_of_Est_cum[Index_prodv1[2]-1]+j2,i);
        }
        for (int k1 = 1; k1 <= Number_of_Est[Index_prodv2[1]]; k1++) {
          for (int k2 = 1; k2 <= Number_of_Est[Index_prodv2[2]]; k2++) {
            for (int i=1; i<=3;i++) {
              Index_r2[i] = J_all (Number_of_Est_cum[Index_prodv2[1]-1]+k1,i);
              Index_r2[i+3] = J_all (Number_of_Est_cum[Index_prodv2[2]-1]+k2,i);
              Index_r12[i] = Index_r1[i];
              Index_r12[i+3] = Index_r1[i+3];
              Index_r12[i+6] = Index_r2[i];
              Index_r12[i+9] = Index_r2[i+3];
            }
            E_1 = Calculate_E_Sr_cpp(Index_S1, Index_r1, prod1plusepssq, gam, K, m, s_squared, v, n);
            E_2 = Calculate_E_Sr_cpp(Index_S2, Index_r2, prod1plusepssq, gam, K, m, s_squared, v, n);
            E_12 = Calculate_E_Sr_cpp(Index_S12, Index_r12, prod1plusepssq, gam, K, m, s_squared, v, n);
            Covariance = Covariance + g(Index_prodv1[1],j1) * g(Index_prodv1[2], j2) * g(Index_prodv2[1],k1) * g(Index_prodv2[2], k2) * (E_12 - E_1 * E_2);
          }
        }
      }
    }

  } else {
    IntegerVector Index_r1 (7);
    IntegerVector Index_r2 (4);
    IntegerVector Index_r12 (10);
    if (Dim1 == 1) {
      int intTemp;
      for (int i = 1; i <= 2; i++) {
        intTemp = Index_S1[i];
        Index_S1[i] = Index_S2[i];
        Index_S2[i] = intTemp;
        intTemp = Index_prodv1[i];
        Index_prodv1[i] = Index_prodv2[i];
        Index_prodv2[i] = intTemp;
      }
    }
    for (int j1 = 1; j1 <= Number_of_Est[Index_prodv1[1]]; j1++) {
      for (int j2 = 1; j2 <= Number_of_Est[Index_prodv1[2]]; j2++) {
        for (int i=1; i<=3;i++) {
          Index_r1[i] = J_all (Number_of_Est_cum[Index_prodv1[1]-1]+j1,i);
          Index_r1[i+3] = J_all (Number_of_Est_cum[Index_prodv1[2]-1]+j2,i);
        }
        for (int k1 = 1; k1 <= Number_of_Est[Index_prodv2[1]]; k1++) {
          for (int i=1; i<=3;i++) {
            Index_r2[i] = J_all (Number_of_Est_cum[Index_prodv2[1]-1]+k1,i);
            Index_r12[i] = Index_r1[i];
            Index_r12[i+3] = Index_r1[i+3];
            Index_r12[i+6] = Index_r2[i];
          }
          E_1 = Calculate_E_Sr_cpp(Index_S1, Index_r1, prod1plusepssq, gam, K, m, s_squared, v, n);
          E_2 = Calculate_E_Sr_cpp(Index_S2, Index_r2, prod1plusepssq, gam, K, m, s_squared, v, n);
          E_12 = Calculate_E_Sr_cpp(Index_S12, Index_r12, prod1plusepssq, gam, K, m, s_squared, v, n);
          Covariance = Covariance + g(Index_prodv1[1],j1) * g(Index_prodv1[2], j2) * g(Index_prodv2[1], k1) * (E_12 - E_1 * E_2);
        }
      }
    }
  }
  //Rcout << E_1 << " " << E_2 << " " << E_12 << "      ";
  //Rcout << Covariance << "\n";
  return(Covariance);
}



// //' @export
// [[Rcpp::export]]
double Cov_S_hat_infty1_S_hat_infty2_cpp(IntegerVector J_s, IntegerVector K_s, NumericVector g_s, NumericVector h_s, IntegerVector Index_S_hat_infty1_s, IntegerVector Index_S_hat_infty2_s, NumericVector epsilon_s, NumericVector m_s, NumericVector s_squared_s, NumericVector v_s, int n, IntegerVector Number_of_Est_s) {
  int i;
  int j;
  int k;
  int l;
  int mu;
  int nu;
  double dblTemp;
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
  //Rcout << "J_all\n";
  //for (i = 1; i <= Number_of_Est_cum[n]; i++) {
  //  for (k = 1; k <= 3; k++) {
  //    Rcout << J_all(i,k) << " ";
  //  }
  //  Rcout << "\n";
  //}


  NumericMatrix g (n+1, Max_Number_of_Est + 1); // g(i,j) entspricht g[[i]][j]
  for (i = 1; i <= n; i++) {
    for (j = Number_of_Est_cum[i-1]; j < Number_of_Est_cum[i]; j++) {
      g(i, j+1-Number_of_Est_cum[i-1]) = g_s[j];
    }
  }
  //Rcout << "g\n";
  //for (i=1; i<=n; i++) {
  //  for (k=1; k<=Max_Number_of_Est; k++) {
  //    Rcout << g(i,k) << " ";
  //  }
  //  Rcout << "\n";
  //}


  NumericMatrix h (n + 1, n + 1); // h(k,i) entspricht h[[k]][i]
  l=0;
  for (k = 1; k <= n; k++) {
    for (i = 1; i < n-k+2; i++) {
      h(k, i) = h_s[l];
      l++;
    }
  }

  //Rcout << "h\n";
  //for (i=1; i<=n; i++) {
  //  for (k=1; k<=n; k++) {
  //    Rcout << h(i,k) << " ";
  //  }
  //  Rcout << "\n";
  //}

  IntegerVector Index_S_hat_infty1 (3);
  for (i=1; i<=2; i++) {
    Index_S_hat_infty1(i)=Index_S_hat_infty1_s(i-1);
  }
  IntegerVector Index_S_hat_infty2 (3);
  for (i=1; i<=2; i++) {
    Index_S_hat_infty2(i)=Index_S_hat_infty2_s(i-1);
  }


  int Count_K = K_s.size();
  IntegerVector K (Count_K + 1);
  for (k=1; k <= Count_K; k++) {
    K[k] = K_s[k-1];
  }
  //Rcout << "K\n";
  //for (i=1; i<=Count_K; i++) {
  //  Rcout << K(i) << " ";
  //}
  //Rcout << "\n";


  NumericVector s_squared (n+1);
  for (k=1; k<=n; k++) {
    s_squared[k] = s_squared_s[k-1];
  }
  NumericVector v (n+1);
  for (i=1; i<=n; i++) {
    v[i] = v_s[i-1];
  }
  NumericVector epsilon (n);
  for (i=1; i<=n-1; i++) {
    epsilon[i] = epsilon_s[i-1];
  }
  NumericVector m (n+1);
  for (k=1; k<=n; k++) {
    m[k] = m_s[k-1];
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

  //Rcout << "gam\n";
  //for (i = 1; i <= n; i++) {
  //  for (k = 1; k <= n; k++) {
  //    Rcout << gam(i,k) << " ";
  //  }
  //  Rcout << "\n";
  //}

  // calculate squared epsilons
  NumericVector epsilon_sqared (n);
  for (i=1; i<=n-1; i++) {
    epsilon_sqared[i] = epsilon[i] * epsilon[i];
  }
  //Rcout << "epsilon_sqared\n";
  //for (i = 1; i <= n-1; i++) {
  //  Rcout << epsilon_sqared[i] << " ";
  //}
  //Rcout << "\n";

  // calculate \prod_{\nu=i}^{n-1} (1+\varepsilon_\nu^2)
  // set to 1 if i=n
  NumericVector prod1plusepssq (n+1);
  prod1plusepssq (n) = 1;
  for (i=n-1; i>=1; i--) {
    prod1plusepssq[i] = prod1plusepssq[i+1] * (1+epsilon_sqared[i]);
  }
  //Rcout << "prod1plusepssq\n";
  //for (i = 1; i <= n; i++) {
  //  Rcout << prod1plusepssq[i] << " ";
  //}
  //Rcout << "\n";



  i = Index_S_hat_infty1[1];
  k = Index_S_hat_infty1[2];
  j = Index_S_hat_infty2[1];
  l = Index_S_hat_infty2[2];

  double Covariance = 0;

  IntegerVector Index_S1 (3);
  IntegerVector Index_S2 (3);
  IntegerVector Index_prodv1 (3);
  IntegerVector Index_prodv2 (3);

  for (mu = 1; mu <=n-k+1; mu++) {
    for (nu = 1; nu <= n-l+1; nu++) {
      dblTemp = 0;
      Index_S1[1] = mu;
      Index_S1[2] = k;
      Index_prodv1[1] = mu;
      Index_prodv1[2] = i;
      Index_S2[1] = nu;
      Index_S2[2] = l;
      Index_prodv2[1] = nu;
      Index_prodv2[2] = j;
      dblTemp = dblTemp + 4 * Cov_S_by_vhat_S_by_vhat_cpp(J_all, g, Index_S1, Index_prodv1, 1, Index_S2, Index_prodv2, 1, prod1plusepssq, gam, K, m, s_squared, v, n, Number_of_Est, Number_of_Est_cum);
      //Rcout << dblTemp << " ";
      dblTemp = dblTemp - 2 * v[i] * Cov_S_by_vhat_S_by_vhat_cpp(J_all, g, Index_S1, Index_prodv1, 2, Index_S2, Index_prodv2, 1, prod1plusepssq, gam, K, m, s_squared, v, n, Number_of_Est, Number_of_Est_cum);
      //Rcout << dblTemp << " ";
      dblTemp = dblTemp - 2 * v[j] * Cov_S_by_vhat_S_by_vhat_cpp(J_all, g, Index_S1, Index_prodv1, 1, Index_S2, Index_prodv2, 2, prod1plusepssq, gam, K, m, s_squared, v, n, Number_of_Est, Number_of_Est_cum);
      //Rcout << dblTemp << " ";
      dblTemp = dblTemp + v[i] * v[j] * Cov_S_by_vhat_S_by_vhat_cpp(J_all, g, Index_S1, Index_prodv1, 2, Index_S2, Index_prodv2, 2, prod1plusepssq, gam, K, m, s_squared, v, n, Number_of_Est, Number_of_Est_cum);
      //Rcout << dblTemp << " ";

      //Rcout << dblTemp;
      Covariance = Covariance + h(k,mu) * h(l,nu) * dblTemp;

      //Rcout << Covariance << "\n";
    }
  }
  //Rcout << versuch(gam);
  //NumericMatrix test (n+1, n+1);
  //test =  versuch(gam);
  //Rcout << "test\n";
  //for (i = 1; i <= n; i++) {
  //  for (k = 1; k <= n; k++) {
  //    Rcout << test(i,k) << " ";
  //  }
  //  Rcout << "\n";
  //}

  //Rcout << g;
  Covariance = v[i] * v[j] * Covariance;
  return Covariance;
}





