

if (Dim1 == 1 && Dim2 == 1) {
  for (int j1 = 1; j1 <= Number_of_Est[Index_prodv1[1]]; j1++) {
    Index_r1 <- matrix(J[[Index_prodv1[1]]][j1,],ncol = 3)
    for (int k1 = 1; k1 <= Number_of_Est[Index_prodv2[1]]; k1++) {
      Index_r2 <- matrix(J[[Index_prodv2[1]]][k1,],ncol = 3)
      Index_r12 <- rbind(Index_r1, Index_r2)
      E_1 <- Calculate_E_Sr_cpp(Index_S1, Index_r1, prod1plusepssq, gam, K, m, s_squared, v)
      E_2 <- Calculate_E_Sr_cpp(Index_S2, Index_r2, prod1plusepssq, gam, K, m, s_squared, v)
      E_12 <- Calculate_E_Sr_cpp(Index_S12, Index_r12, prod1plusepssq, gam, K, m, s_squared, v)
      Covariance <- Covariance + g[[Index_prodv1[1]]][j1] * g[[Index_prodv2[1]]][k1] * (E_12 - E_1 * E_2)
    }
  }
} else if (length(Index_prodv1) == 2 && length(Index_prodv2) == 2) {
  for (j1 in 1:nrow(J[[Index_prodv1[1]]])) {
    for (j2 in 1:nrow(J[[Index_prodv1[2]]])) {
      Index_r1 <- matrix(c(J[[Index_prodv1[1]]][j1,],J[[Index_prodv1[2]]][j2,]),ncol = 3,byrow = T)
      for (k1 in 1:nrow(J[[Index_prodv2[1]]])) {
        for (k2 in 1:nrow(J[[Index_prodv2[2]]])) {
          Index_r2 <- matrix(c(J[[Index_prodv2[1]]][k1,],J[[Index_prodv2[2]]][k2,]), ncol = 3, byrow = T)
          Index_r12 <- rbind(Index_r1, Index_r2)
          E_1 <- Calculate_E_Sr_cpp(Index_S1, Index_r1, prod1plusepssq, gam, K, m, s_squared, v)
          E_2 <- Calculate_E_Sr_cpp(Index_S2, Index_r2, prod1plusepssq, gam, K, m, s_squared, v)
          E_12 <- Calculate_E_Sr_cpp(Index_S12, Index_r12, prod1plusepssq, gam, K, m, s_squared, v)
          Covariance <- Covariance + g[[Index_prodv1[1]]][j1] * g[[Index_prodv1[2]]][j2] * g[[Index_prodv2[1]]][k1] * g[[Index_prodv2[2]]][k2] * (E_12 - E_1 * E_2)
        }
      }
    }
  }

} else {
  if (length(Index_prodv1) == 1) {
    Temp <- Index_S1
    Index_S1 <- Index_S2
    Index_S2 <- Temp
    Temp <- Index_prodv1
    Index_prodv1 <- Index_prodv2
    Index_prodv2 <- Temp
  }
  for (j1 in 1:nrow(J[[Index_prodv1[1]]])) {
    for (j2 in 1:nrow(J[[Index_prodv1[2]]])) {
      Index_r1 <- matrix(c(J[[Index_prodv1[1]]][j1,],J[[Index_prodv1[2]]][j2,]),ncol = 3,byrow = T)
      for (k1 in 1:nrow(J[[Index_prodv2[1]]])) {
        Index_r2 <- matrix(J[[Index_prodv2[1]]][k1,], ncol = 3, byrow = T)
        Index_r12 <- rbind(Index_r1, Index_r2)
        E_1 <- Calculate_E_Sr_cpp(Index_S1, Index_r1, prod1plusepssq, gam, K, m, s_squared, v)
        E_2 <- Calculate_E_Sr_cpp(Index_S2, Index_r2, prod1plusepssq, gam, K, m, s_squared, v)
        E_12 <- Calculate_E_Sr_cpp(Index_S12, Index_r12, prod1plusepssq, gam, K, m, s_squared, v)
        Covariance <- Covariance + g[[Index_prodv1[1]]][j1] * g[[Index_prodv1[2]]][j2] * g[[Index_prodv2[1]]][k1] * (E_12 - E_1 * E_2)
      }
    }
  }
}
return(Covariance)
