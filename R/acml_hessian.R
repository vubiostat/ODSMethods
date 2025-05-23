# Lucy's Analytical Hessian

transform_output <- function(InitVals, x = x, y = y, z = z){
  nx <- ncol(x)
  beta = InitVals[1:nx] # note here are with beta_0, need to change to accommodate if no beta_0
  log_sigma_vc = InitVals[(nx+1):(nx+ncol(z))]
  log_inv_rho_vc = InitVals[nx+3]
  log_sigma_e = InitVals[length(InitVals)]
  sigma_vc = exp(log_sigma_vc)
  sigma_e = exp(log_sigma_e)
  rho_vc = tanh(log_inv_rho_vc)
  list(beta = beta, sigma_vc = sigma_vc, rho_vc = rho_vc, sigma_e = sigma_e)
}

vi_calc = function(zi, sigma_vc, rho_vc, sigma_e){
  SDMat_RE  = diag(sigma_vc)
  ncolzi    = ncol(zi) ## make sure this equals length(sigma.vc)
  nrowzi    = nrow(zi)
  nERRsd    = length(sigma_e)
  b         = matrix(0,ncolzi,ncolzi)
  b[lower.tri(b, diag=FALSE)] = rho_vc
  CorMat_RE = t(b)+b+diag(rep(1,ncolzi))
  CovMat_RE = SDMat_RE %*% CorMat_RE %*% SDMat_RE
  zi %*% CovMat_RE %*% t(zi) + diag(rep(sigma_e^2, each=nrowzi/nERRsd))
}

# Not exported, needs interface work and validation
acml_hessian <- function(
    estimate, y, x, z, id, cutpoints, w_function,
    Weights  = rep(1,length(y)),
    SampProb = rep(1, length(cutpoints)+1))
{
  output_org = transform_output(estimate, x = x, y = y, z = z)
  beta = output_org$beta
  log_sigma_vc = estimate[(ncol(x)+2):(ncol(x)+ncol(z)+2)]
  log_inv_rho_vc = estimate[(ncol(x)+ncol(z)+3):(ncol(x)+ncol(z) + choose(ncol(z)+1, 2)+2)]
  log_sigma_e = estimate[(ncol(x)+ncol(z) + choose(ncol(z)+1, 2)+3):length(estimate)]
  sigma_vc = output_org$sigma_vc
  sigma_e = output_org$sigma_e
  rho_vc = output_org$rho_vc

  Info_mat = Reduce("+",lapply(unique(id),function(i){
    yi = y[id == i]
    Weight_i = unique(Weights[id == i])

    xi = as.matrix(x[id == i,])
    zi = as.matrix(z[id == i,])
    resid   = yi - xi %*% beta
    vi      = vi_calc(zi, sigma_vc, rho_vc, sigma_e)
    inv_vi  = solve(vi)
    mu      = xi %*% beta
    if (w_function == "intercept"){
      wi      = matrix((solve(t(zi[,1:2])%*% zi[,1:2]) %*% t(zi[,1:2]))[1,], nrow = 1) # Intercept
    } else if (w_function == "slope"){
      wi      = matrix((solve(t(zi[,1:2])%*% zi[,1:2]) %*% t(zi[,1:2]))[2,], nrow = 1) # Slope
    }

    mu_q    = wi %*% mu
    Sigma_q = wi %*% vi %*% t(wi)

    b0         = matrix(0, length(sigma_vc), length(sigma_vc))
    b0[lower.tri(b0, diag=FALSE)] = rho_vc
    tbb0       = t(b0) + b0
    CorMat.RE  = tbb0 +  diag(rep(1,length(sigma_vc)))

    D             = diag(sigma_vc) %*% CorMat.RE %*% diag(sigma_vc)
    m             = matrix(0,length(sigma_vc), length(sigma_vc))

    ## derivatives w.r.t variance components SDs
    dVMat_dlog_sigma_vc = NULL

    for (mmm in 1:length(sigma_vc)){
      m1               = m
      m1[mmm,]         = m1[,mmm] = 1
      m1[mmm,mmm]      = 2
      tmp              = sigma_vc[mmm]+ 1*(sigma_vc[mmm]==0) ## to prevent unlikely division by 0
      dDMat_dsigma_vc  = m1 * D / tmp
      dVMat_dlog_sigma_vc[[mmm]]  = zi %*% dDMat_dsigma_vc %*% t(zi) * tmp
    }

    ## derivatives w.r.t variance components rhos
    b1         = matrix(0,length(sigma_vc),length(sigma_vc))
    b1[lower.tri(b1, diag=FALSE)] =  c(1:length(rho_vc))
    tbb1       = t(b1) + b1
    dVMat_dtrans_rho_vc = NULL

    for (mmm in 1:length(rho_vc)){
      tmp                   = which(tbb1==mmm, arr.ind=TRUE)
      m2                    = m
      m2[tmp[1,1],tmp[1,2]] = 1
      m2[tmp[2,1],tmp[2,2]] = 1
      tmp                   = rho_vc[mmm] + 1*(rho_vc[mmm]==0) ## to prevent division by 0
      dDMat_drho_vc         = m2 * D / tmp
      dVMat_dtrans_rho_vc[[mmm]]         = zi %*% dDMat_drho_vc %*% t(zi) * ((1-tmp^2)/2)
    }

    ## derivatives w.r.t error sds
    dVMat_dlog_sigma_e = NULL

    for (mmm in 1:length(sigma_e)){
      dsigma_e_vec       = 2*sigma_e
      dsigma_e_vec[-mmm] = 0
      dVMat_dlog_sigma_e[[mmm]]    = diag(rep(dsigma_e_vec * sigma_e[mmm], each=length(yi)/length(sigma_e)))
    }

    ## second derivatives w.r.t error sds
    d2VMat_dlog_sg_e_dlog_sg_e = NULL
    for (mmm in 1:length(sigma_e)){
      dsigma_e_vec       = 4
      d2VMat_dlog_sg_e_dlog_sg_e[[mmm]]    = diag(rep(dsigma_e_vec*sigma_e[mmm]^2, each=length(yi)/length(sigma_e)))
    }

    d2VMat_dlog_sg_e_dalpha = NULL
    for (mmm in 1:length(sigma_e)){
      d2VMat_dlog_sg_e_dalpha[[mmm]]    = diag(rep(0, each=length(yi)/length(sigma_e)))
    }

    ## second derivatives w.r.t variance components SDs
    d2VMat_dlog_sg_vc_dlog_sg_vc = NULL
    for (mmm in 1:length(sigma_vc)){
      d2VMat_dlog_sg_vc_dlog_sg_vc_tmp = NULL
      for (nnn in 1:length(sigma_vc)){
        if (mmm == nnn){
          m1               = m
          m1[mmm,]         = m1[,mmm] = 1
          m1[mmm,mmm]      = 4
          d2DMat_dsigma_vc2  = m1 * D
          d2VMat_dlog_sg_vc_dlog_sg_vc_tmp[[nnn]]  = zi %*% (d2DMat_dsigma_vc2) %*% t(zi)

        } else {
          m1               = m
          m1[mmm,nnn]         = m1[nnn,mmm] = 1
          tmp1              = sigma_vc[mmm]+ 1*(sigma_vc[mmm]==0) ## to prevent unlikely division by 0
          tmp2              = sigma_vc[nnn]+ 1*(sigma_vc[nnn]==0) ## to prevent unlikely division by 0
          d2DMat_dalpha_dalpha  = m1 * D
          d2VMat_dlog_sg_vc_dlog_sg_vc_tmp[[nnn]]  = zi %*% d2DMat_dalpha_dalpha %*% t(zi)

        }
      }
      d2VMat_dlog_sg_vc_dlog_sg_vc[[mmm]] = d2VMat_dlog_sg_vc_dlog_sg_vc_tmp
    }

    ## second derivatives w.r.t variance components rhos
    b1         = matrix(0,length(sigma_vc),length(sigma_vc))
    b1[lower.tri(b1, diag=FALSE)] =  c(1:length(rho_vc))
    tbb1       = t(b1) + b1
    d2VMat_dtrans_rho_vc_dtrans_rho_vc = NULL
    for (mmm in 1:length(rho_vc)){
      d2VMat_dtrans_rho_vc_dtrans_rho_vc_tmp = NULL
      for (nnn in 1:length(rho_vc)){
        if (mmm == nnn){
          tmp                   = which(tbb1==mmm, arr.ind=TRUE)
          m2                    = m
          m2[tmp[1,1],tmp[1,2]] = 1
          m2[tmp[2,1],tmp[2,2]] = 1
          tmp1                   = rho_vc[mmm] + 1*(rho_vc[mmm]==0) ## to prevent division by 0
          # tmp2                   = rho_vc[nnn] + 1*(rho_vc[nnn]==0) ## to prevent division by 0
          dDMat_drho_vc         = m2 * D / tmp1
          d2VMat_dtrans_rho_vc_dtrans_rho_vc_tmp[[nnn]]         = zi %*% dDMat_drho_vc %*% t(zi) * ((tmp1^3 - tmp1)/2)
        } else {
          tmp                   = which(tbb1==mmm, arr.ind=TRUE)
          m2                    = m
          m2[tmp[1,1],tmp[1,2]] = 0
          m2[tmp[2,1],tmp[2,2]] = 0
          tmp1                   = rho_vc[mmm] + 1*(rho_vc[mmm]==0) ## to prevent division by 0
          tmp2                   = rho_vc[nnn] + 1*(rho_vc[nnn]==0) ## to prevent division by 0
          dDMat_drho_vc         = m2 * D / tmp1
          d2VMat_dtrans_rho_vc_dtrans_rho_vc_tmp[[nnn]]         = zi %*% dDMat_drho_vc %*% t(zi)
        }
      }
      d2VMat_dtrans_rho_vc_dtrans_rho_vc[[mmm]] = d2VMat_dtrans_rho_vc_dtrans_rho_vc_tmp
    }


    ## second derivatives w.r.t variance components SDs
    b1         = matrix(0,length(sigma_vc),length(sigma_vc))
    b1[lower.tri(b1, diag=FALSE)] =  c(1:length(rho_vc))
    tbb1       = t(b1) + b1
    d2VMat_dlog_sg_vc_dtrans_rho_vc = NULL
    for (mmm in 1:length(sigma_vc)){
      d2VMat_dlog_sg_vc_dtrans_rho_vc_tmp = NULL
      for (nnn in 1:length(rho_vc)){
        rho_loc                       = which(tbb1==nnn, arr.ind=TRUE)
        m2                            = m
        m2[rho_loc[1,1],rho_loc[1,2]] = 1
        m2[rho_loc[2,1],rho_loc[2,2]] = 1
        tmp1                          = sigma_vc[mmm]+ 1*(sigma_vc[mmm]==0) ## to prevent unlikely division by 0
        tmp2                          = rho_vc[nnn] + 1*(rho_vc[nnn]==0) ## to prevent division by 0
        if (mmm %in% rho_loc) {
          d2DMat_dalpha_dalpha = m2 * D / tmp1 / tmp2
        } else {
          d2DMat_dalpha_dalpha = diag(0, nrow=nrow(m2),ncol = ncol(m2))
        }
        d2VMat_dlog_sg_vc_dtrans_rho_vc_tmp[[nnn]]  = zi %*% d2DMat_dalpha_dalpha %*% t(zi) * tmp1 * (1-tmp2^2)/2
        # print(d2DMat_dalpha_dalpha)
      }
      d2VMat_dlog_sg_vc_dtrans_rho_vc[[mmm]] = d2VMat_dlog_sg_vc_dtrans_rho_vc_tmp
    }

    tot_length = length(sigma_vc)+length(rho_vc)+length(sigma_e)
    dVMat_all_list <- c(dVMat_dlog_sigma_vc, dVMat_dtrans_rho_vc, dVMat_dlog_sigma_e)

    # if (sum(dim(Sigma_q) != c(1,1)) == 0){
    Sigma_q = Sigma_q[1,1]
    resid = yi - xi %*% beta
    ACi = ACi1q_new(cutpoints, SampProb, mu_q, Sigma_q)
    diff_pi = SampProb[1:(length(SampProb)-1)] - SampProb[2:(length(SampProb))]
    Phi = dnorm((cutpoints - c(mu_q))*(Sigma_q^(-0.5)),0,1)
    d_Phi = d_std_norm_density((cutpoints - c(mu_q))*(Sigma_q^(-0.5)))
    # test_mult = c(log_sigma_vc,(-rho_vc)*(1-rho_vc^2)/2,log_sigma_e)

    if (sum(SampProb != 1) == 0) {
      I1_i = t(xi) %*% solve(vi) %*% xi * Weight_i
    } else {
      I1_i = t(xi) %*% solve(vi) %*% xi - ACi^(-2)*Sigma_q^(-1)*t(wi %*% xi)%*%(wi %*% xi)*(sum(diff_pi*Phi))^2 + ACi^(-1)*Sigma_q^(-1)*t(wi %*% xi)%*%(wi %*% xi)*sum(diff_pi*d_Phi)
    }

    I4_i <- matrix(NA, nrow = tot_length,ncol =tot_length)
    if (sum(SampProb != 1) == 0) {
      for (ii in 1:tot_length){
        for (jj in 1:tot_length){
          I4_i[ii,jj] = 0.5*sum(diag(inv_vi %*% dVMat_all_list[[ii]] %*% inv_vi %*% dVMat_all_list[[jj]] * Weight_i))
        }
      }
    } else {
      l_sigma_vc = length(sigma_vc)
      l_rho_vc = length(rho_vc)
      l_sigma_e = length(sigma_e)
      loc_sigma_vc = 1:length(sigma_vc)
      loc_rho_vc = (length(sigma_vc)+1):(length(sigma_vc)+length(rho_vc))
      loc_sigma_e = (length(sigma_vc)+length(rho_vc)+1):tot_length

      # dli_dalpha_all_list = c(unlist(dli_dsigma_vc), unlist(dli_drho_vc), unlist(dli_dsigma_e))

      for (ii in 1:tot_length){
        for (jj in 1:tot_length){
          if ((ii %in% loc_sigma_e) & (jj %in% loc_sigma_e) & ii == jj){
            d2VMat_dalpha_dalpha = d2VMat_dlog_sg_e_dlog_sg_e[[1]]
          } else if ((ii %in% loc_sigma_e) | (jj %in% loc_sigma_e)){
            d2VMat_dalpha_dalpha = d2VMat_dlog_sg_e_dalpha[[1]]
          } else if ((ii %in% loc_sigma_vc) & (jj %in% loc_sigma_vc)) {
            d2VMat_dalpha_dalpha = d2VMat_dlog_sg_vc_dlog_sg_vc[[ii]][[jj]]
          } else if ((ii %in% loc_rho_vc) & (jj %in% loc_rho_vc)) {
            d2VMat_dalpha_dalpha = d2VMat_dtrans_rho_vc_dtrans_rho_vc[[ii-l_sigma_vc]][[jj-l_sigma_vc]]
          } else if ((ii %in% loc_sigma_vc) & (jj %in% loc_rho_vc)) {
            d2VMat_dalpha_dalpha = d2VMat_dlog_sg_vc_dtrans_rho_vc[[ii]][[jj-l_sigma_vc]]
          } else if ((jj %in% loc_sigma_vc) & (ii %in% loc_rho_vc)) {
            d2VMat_dalpha_dalpha = d2VMat_dlog_sg_vc_dtrans_rho_vc[[jj]][[ii-l_sigma_vc]]
          }

          I4_i[ii,jj] =  0.5*sum(diag(inv_vi%*%d2VMat_dalpha_dalpha-inv_vi%*%dVMat_all_list[[ii]] %*% inv_vi %*% dVMat_all_list[[jj]])) - 0.5*t(resid)%*%inv_vi%*%(d2VMat_dalpha_dalpha - 2*dVMat_all_list[[ii]] %*% inv_vi %*% dVMat_all_list[[jj]])%*%inv_vi%*%resid - (1/4)*ACi^(-2)*Sigma_q^(-3)* t(wi%*%dVMat_all_list[[ii]]%*%t(wi)) %*% (wi%*%dVMat_all_list[[jj]]%*%t(wi)) * (sum((cutpoints - c(mu_q))*diff_pi*Phi))^2  + (3/4)*ACi^(-1)*Sigma_q^(-5/2) *t(wi%*%dVMat_all_list[[ii]]%*%t(wi)) %*% (wi%*%dVMat_all_list[[jj]]%*%t(wi))  * sum((cutpoints - c(mu_q))*diff_pi*Phi) - (1/2)*ACi^(-1)*Sigma_q^(-3/2) * (wi%*%d2VMat_dalpha_dalpha%*%t(wi)) * sum((cutpoints - c(mu_q))*diff_pi*Phi) + (1/4)*ACi^(-1)*Sigma_q^(-3) * t(wi%*%dVMat_all_list[[ii]]%*%t(wi)) %*% (wi%*%dVMat_all_list[[jj]]%*%t(wi)) * sum((cutpoints - c(mu_q))^2*diff_pi*d_Phi)
        }
      }
    }
    if (sum(SampProb != 1) == 0) {
      I3_i <- matrix(0, nrow = length(beta), ncol = tot_length)
    } else {
      I3_i <- matrix(NA, nrow = length(beta),ncol = tot_length)
      for (iii in 1:tot_length){
        I3_i[,iii] = t(xi)%*%inv_vi%*%dVMat_all_list[[iii]]%*%inv_vi%*%resid + t(- (1/2)*ACi^(-2)*Sigma_q^(-2)*(wi%*%dVMat_all_list[[iii]]%*%t(wi))%*%(wi%*%xi)*sum((cutpoints - c(mu_q))*diff_pi*Phi)*sum(diff_pi*Phi) + (1/2)*ACi^(-1)*Sigma_q^(-3/2)*(wi%*%dVMat_all_list[[iii]]%*%t(wi))%*%(wi%*%xi)*sum(diff_pi*Phi) + (1/2)*ACi^(-1)*Sigma_q^(-2)*(wi%*%dVMat_all_list[[iii]]%*%t(wi))%*%(wi%*%xi)*sum((cutpoints - c(mu_q))*diff_pi*d_Phi))
      }
    }
    I_i = matrix(NA, nrow = length(estimate),ncol = length(estimate))
    I_i[1:length(beta),1:length(beta)] = I1_i
    I_i[1:length(beta),(length(beta)+1):length(estimate)] = I3_i
    I_i[(length(beta)+1):length(estimate),1:length(beta)] = t(I3_i)
    I_i[(length(beta)+1):length(estimate),(length(beta)+1):length(estimate)] = I4_i
    I_i
  }))
  -Info_mat
}
