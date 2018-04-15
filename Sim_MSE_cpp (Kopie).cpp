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

// [[Rcpp::export]]
NumericVector Sim_MSE_bullet_GLR_cpp(NumericVector S_sim, NumericVector v_sim, int n, long NumberOfSimulations, NumericVector K_s, NumericVector J_s, IntegerVector Number_of_Est_s, NumericVector h_s, NumericVector g_s, NumericVector m_s, NumericVector s_squared_s) {

  Rcpp::checkUserInterrupt();

  int i;
  int j;
  int k;
  int l;
  long SimNr;
  long Counter_v = 0;
  long Counter_S = 0;
  NumericVector v(n+1);
  NumericMatrix S(n+1,n+1);

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

  NumericMatrix h (n + 1, n + 1); // h(k,i) entspricht h[[k]][i]
  l=0;
  for (k = 1; k <= n; k++) {
    for (i = 1; i < n-k+2; i++) {
      h(k, i) = h_s[l];
      l++;
    }
  }

  NumericVector m (n+1);
  for (k=1; k<=n; k++) {
    m[k] = m_s[k-1];
  }


  NumericVector s_squared (n+1);
  for (k=1; k<=n; k++) {
    s_squared[k] = s_squared_s[k-1];
  }


  // calculate gammas
  // set gam[i,l] <- 1 if $l\not\in K$

  NumericMatrix gam (n + 1, n + 1);
  for (i = 1; i <= n; i++) {
    for (k = 1; k <= n; k++) {
      gam(i,k) = 1;
    }
  }



  NumericMatrix r_hat (n + 1, Max_Number_of_Est + 1);
   double Temp;
   NumericVector v_hat (n+1);
   NumericVector m_hat(n+1);
   NumericVector Reserve(n+1);
   NumericVector Reserve_hat(n+1);
   NumericVector SquaredError (n+1);

   for (i = 0; i <= n; i++) {
     SquaredError[i] = 0;
   }


  for (SimNr = 1; SimNr <= NumberOfSimulations; SimNr++) {
    for (i=1; i<=n; i++) {
      v[i]= v_sim[Counter_v];
      Counter_v++;
    }
    for (k = 1; k <=n; k++) {
      for (i = 1; i <= n; i++) {
        S(i,k) = S_sim[Counter_S];
        Counter_S++;
      }
    }
    for (i = 1; i <= n; i++) {
      for (k = 1; k <= Count_K; k++) {
        gam(i,K[k]) = 1+s_squared[K[k]]/(v[i]*m[K[k]]*m[K[k]]);
      }
    }

    // Calculate estimators r_hat
    for (i = 1; i <= n; i++) {
      for (int nu = 1; nu <= Number_of_Est[i]; nu++) {
        j = J_all(Number_of_Est_cum[i-1]+nu, 2);
        k = J_all(Number_of_Est_cum[i-1]+nu, 3);
        if (j == 0) {
          r_hat(i,nu) = 1/v[i];
        } else {
          r_hat(i,nu) = 1/v[j] * S(j,k) / S(i,k) / gam(i,k);
        }
      }
    }

    for (i = 1; i <= n; i++) {
      Temp = 0;
      for (j = 1; j <= Number_of_Est[i]; j++) {
        Temp += g(i,j) * r_hat(i,j);
      }
      v_hat[i] = 1 / Temp;
    }

    for (k = 1; k <= n; k++) {
      m_hat[k] = 0;
      for (i = 1; i <= n-k+1; i++) {
        m_hat[k] += h(k,i) * S(i,k) / v_hat[i];
      }
    }

    Reserve[0] = 0;
    Reserve_hat[0] = 0;
    Reserve[n] = 0;
    Reserve_hat[n] = 0;
    for (i = 2; i <= n; i++) {
      Reserve[i-1] = 0;
      Reserve_hat[i-1] = 0;
      for (k = n-i+2; k <= n; k++) {
        Reserve[i-1] += S(i,k);
        Reserve_hat[i-1] += v_hat[i] * m_hat[k];
      }
      Reserve[n] += Reserve[i-1];
      Reserve_hat[n] += Reserve_hat[i-1];
    }

    for (i = 0; i <= n; i++) {
      SquaredError[i] += (Reserve[i] - Reserve_hat[i]) * (Reserve[i] - Reserve_hat[i]);
    }
    //Rcout << SimNr << " ";
  }

  //Rcout << "\nr_hat\n";
  //for (i = 1; i <=n; i++) {
  //  Rcout << "\n";
  ////  for (k=1; k<=Number_of_Est[i]; k++) {
  //  Rcout << r_hat(i,k) << " ";
  //  }
  //}
  //Rcout << " ";

  //Rcout << "\ngam\n";
  //for (i = 1; i <=n; i++) {
  //   Rcout << "\n";
  //   for (k=1; k<=n; k++) {
  //     Rcout << gam(i,k) << " ";
  //   }
  //}
  //Rcout << " ";

  //Rcout << "\n\n v and v_hat\n";
  //for (i = 1; i <=n; i++) {
  //   Rcout << "\n";
  //   Rcout << v[i] << " " << v_hat[i];
  //}
  //Rcout << " ";


  for (i = 0; i <= n; i++) {
    SquaredError[i] = SquaredError[i] / NumberOfSimulations;
  }




  return SquaredError;
}






// [[Rcpp::export]]
NumericVector Sim_MSE_bullet_LR_cpp(NumericVector S_sim, NumericVector v_sim, int n, long NumberOfSimulations, NumericVector W_s) {

  Rcpp::checkUserInterrupt();


  int i;
  int j;
  int k;
  int l;
  long SimNr;
  long Counter_v = 0;
  long Counter_S = 0;
  NumericVector v(n+1);
  NumericMatrix S(n+1,n+1);



  NumericMatrix W (n + 1, n + 1);
  l=0;
  for (k = 1; k <= n; k++) {
    for (i = 1; i <= n; i++) {
      W(i, k) = W_s[l];
      l++;
    }
  }
  double Temp;
  for (k=1; k<= n; k++) {
    Temp = 0;
    for (i=1; i <= n-k+1; i++) {
      Temp += W(i,k);
    }
    for (i=1; i <= n-k+1; i++) {
      W(i,k) = W(i,k) / Temp;
    }
  }



  NumericVector m_hat(n+1);
  NumericVector Reserve(n+1);
  NumericVector Reserve_hat(n+1);
  NumericVector SquaredError (n+1);

  for (i = 0; i <= n; i++) {
    SquaredError[i] = 0;
  }


  for (SimNr = 1; SimNr <= NumberOfSimulations; SimNr++) {
    for (i=1; i<=n; i++) {
      v[i]= v_sim[Counter_v];
      Counter_v++;
    }
    for (k = 1; k <=n; k++) {
      for (i = 1; i <= n; i++) {
        S(i,k) = S_sim[Counter_S];
        Counter_S++;
      }
    }
    for (k = 1; k <= n; k++) {
      m_hat[k] = 0;
      for (i = 1; i <= n-k+1; i++) {
        m_hat[k] += W(i,k) * S(i,k) / v[i];
      }
    }


    Reserve[0] = 0;
    Reserve_hat[0] = 0;
    Reserve[n] = 0;
    Reserve_hat[n] = 0;
    for (i = 2; i <= n; i++) {
      Reserve[i-1] = 0;
      Reserve_hat[i-1] = 0;
      for (k = n-i+2; k <= n; k++) {
        Reserve[i-1] += S(i,k);
        Reserve_hat[i-1] += v[i] * m_hat[k];
      }
      Reserve[n] += Reserve[i-1];
      Reserve_hat[n] += Reserve_hat[i-1];
    }

    for (i = 0; i <= n; i++) {
      SquaredError[i] += (Reserve[i] - Reserve_hat[i]) * (Reserve[i] - Reserve_hat[i]);
    }
    //Rcout << SimNr << " ";
  }

  //Rcout << "\nr_hat\n";
  //for (i = 1; i <=n; i++) {
  //  Rcout << "\n";
  ////  for (k=1; k<=Number_of_Est[i]; k++) {
  //  Rcout << r_hat(i,k) << " ";
  //  }
  //}
  //Rcout << " ";

  //Rcout << "\ngam\n";
  //for (i = 1; i <=n; i++) {
  //   Rcout << "\n";
  //   for (k=1; k<=n; k++) {
  //     Rcout << gam(i,k) << " ";
  //   }
  //}
  //Rcout << " ";

  //Rcout << "\n\n v and v_hat\n";
  //for (i = 1; i <=n; i++) {
  //   Rcout << "\n";
  //   Rcout << v[i] << " " << v_hat[i];
  //}
  //Rcout << " ";


  for (i = 0; i <= n; i++) {
    SquaredError[i] = SquaredError[i] / NumberOfSimulations;
  }




  return SquaredError;
}


// [[Rcpp::export]]
NumericVector Sim_MSE_LR_canonical_cpp(NumericVector S_sim, NumericVector v_sim, int n, long NumberOfSimulations) {

  Rcpp::checkUserInterrupt();


  int i;
  int j;
  int k;
  int l;
  long SimNr;
  long Counter_v = 0;
  long Counter_S = 0;
  NumericVector v(n+1);
  NumericMatrix S(n+1,n+1);



  NumericMatrix W (n + 1, n + 1);
  double Temp;




  NumericVector m_hat(n+1);
  NumericVector Reserve(n+1);
  NumericVector Reserve_hat(n+1);
  NumericVector SquaredError (n+1);

  for (i = 0; i <= n; i++) {
    SquaredError[i] = 0;
  }


  for (SimNr = 1; SimNr <= NumberOfSimulations; SimNr++) {
    for (i=1; i<=n; i++) {
      v[i]= v_sim[Counter_v];
      Counter_v++;
    }

    for (k=1; k<= n; k++) {
      Temp = 0;
      for (i=1; i <= n-k+1; i++) {
        Temp += v[i];
      }
      for (i=1; i <= n-k+1; i++) {
        W(i,k) = v[i] / Temp;
      }
    }

    for (k = 1; k <=n; k++) {
      for (i = 1; i <= n; i++) {
        S(i,k) = S_sim[Counter_S];
        Counter_S++;
      }
    }
    for (k = 1; k <= n; k++) {
      m_hat[k] = 0;
      for (i = 1; i <= n-k+1; i++) {
        m_hat[k] += W(i,k) * S(i,k) / v[i];
      }
    }


    Reserve[0] = 0;
    Reserve_hat[0] = 0;
    Reserve[n] = 0;
    Reserve_hat[n] = 0;
    for (i = 2; i <= n; i++) {
      Reserve[i-1] = 0;
      Reserve_hat[i-1] = 0;
      for (k = n-i+2; k <= n; k++) {
        Reserve[i-1] += S(i,k);
        Reserve_hat[i-1] += v[i] * m_hat[k];
      }
      Reserve[n] += Reserve[i-1];
      Reserve_hat[n] += Reserve_hat[i-1];
    }

    for (i = 0; i <= n; i++) {
      SquaredError[i] += (Reserve[i] - Reserve_hat[i]) * (Reserve[i] - Reserve_hat[i]);
    }
    //Rcout << SimNr << " ";
  }

  //Rcout << "\nr_hat\n";
  //for (i = 1; i <=n; i++) {
  //  Rcout << "\n";
  ////  for (k=1; k<=Number_of_Est[i]; k++) {
  //  Rcout << r_hat(i,k) << " ";
  //  }
  //}
  //Rcout << " ";

  //Rcout << "\ngam\n";
  //for (i = 1; i <=n; i++) {
  //   Rcout << "\n";
  //   for (k=1; k<=n; k++) {
  //     Rcout << gam(i,k) << " ";
  //   }
  //}
  //Rcout << " ";

  //Rcout << "\n\n v and v_hat\n";
  //for (i = 1; i <=n; i++) {
  //   Rcout << "\n";
  //   Rcout << v[i] << " " << v_hat[i];
  //}
  //Rcout << " ";


  for (i = 0; i <= n; i++) {
    SquaredError[i] = SquaredError[i] / NumberOfSimulations;
  }




  return SquaredError;
}
