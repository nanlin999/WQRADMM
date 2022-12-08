#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;

//The soft-thresholding solution of r
arma::vec shrinkcpp(arma::vec u, arma::vec v){
  arma::vec w = (1+sign(u-v))/2%(u-v)-(1+sign(-u-v))/2%(-u-v);
  return w;
}

//The first order derivative of checkloss function 
arma::vec lossweight(arma::vec u, double tau){
  int num = u.size();
  arma::vec loss = arma::zeros(num);
  for(int i = 0; i< num; i++){
    if(u(i) > 0)
      loss(i) = -tau;
    else loss(i) = 1-tau;
  }
  return loss;
}

// The WQR-ADMM algorithm
arma::vec WQRCPP(arma::mat x, arma::vec y, double tau, double rho, double eps, int maxstep, bool intercept, bool warmstart){

  int nsum = x.n_rows;
  if(intercept){
    x.insert_cols(0, arma::ones(nsum));
  }
  int p = x.n_cols;
  arma::vec betaini = arma::zeros(p), zini = arma::zeros(p), uini = arma::zeros(p), vini = arma::zeros(nsum); 
  arma::vec beta = arma::zeros(p), r = arma::zeros(nsum), z = arma::zeros(p), u = arma::zeros(p), v = arma::zeros(nsum); 
  arma::vec xbeta = arma::zeros(nsum); 
  double time = 0, time_warm = 0, time_prep = 0, time_iter = 0;
  arma::vec final = arma::zeros(p+2);
  
  if(warmstart){
    auto start_warm = std::chrono::high_resolution_clock::now();
    zini = inv(x.t()*x)*x.t()*y;
    xbeta = x*zini;
    auto finish_warm = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_warm = finish_warm - start_warm;
    time_warm = elapsed_warm.count();
    time = time + time_warm;
  }
  
  arma::mat tmp = arma::zeros(p, p);
  auto start_prep = std::chrono::high_resolution_clock::now();
  if(nsum>p) 
    tmp = inv(x.t()*x+arma::eye(p,p));
  else tmp = arma::eye(p,p)-x.t()*inv(x*x.t()+arma::eye(nsum,nsum))*x;
  auto finish_prep = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_prep = finish_prep - start_prep;
  time_prep = elapsed_prep.count();
  time = time + time_prep;
  
  int iteration = 0;
  double distance = 1;
  
  while((distance>eps)&&(iteration<maxstep)){
  
    auto start_iter = std::chrono::high_resolution_clock::now();
    //update beta
    beta = zini+uini/rho;
    //update r
    r = shrinkcpp(vini/rho+y-xbeta-0.5*(2*tau-1)/(rho),0.5*arma::ones<arma::vec>(nsum)/(rho));
    //update z
    z = tmp*(x.t()*(y-r+vini/rho)+beta-uini/rho);
    //update u and v
    xbeta = x*z;
    u = uini+rho*(z-beta);
    v = vini+rho*(y-xbeta-r);
    auto finish_iter = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_iter = finish_iter - start_iter;
    time_iter = elapsed_iter.count();
    time = time + time_iter;
    if(accu(abs(betaini))==0)
      distance = 1;
    else distance = accu(abs(beta-betaini))/accu(abs(betaini));
    betaini = beta, zini = z, uini = u, vini = v;
    iteration = iteration + 1;
    
  }
  
  final.subvec(0,p-1) = beta;
  final(p) = time;
  final(p+1) = iteration;
  return final;
  
}

// The adjusted weight
arma::vec WeightCPP(arma::mat x, arma::vec y, arma::vec rep, double tau, arma::vec betahat, double eps, double maxstep, bool intercept, String type){
  
  Environment stats("package:stats");
  Function dnorm = stats["dnorm"];
  Environment base("package:base");
  Function sweep = base["sweep"];
  
  int nsum = x.n_rows;
  if(intercept){
    x.insert_cols(0, arma::ones(nsum));
  }
  int p = x.n_cols, n = rep.size();
  double time_prep1 = 0, time_prep2 = 0, time_iter = 0, time_wei = 0, time = 0;
  
  //calculate D1 and D2
  arma::vec e = arma::zeros(nsum);
  arma::mat xM = arma::zeros(p,nsum), D1 = arma::zeros(p,p), D2 = arma::zeros(p,p);
  auto start_prep1 = std::chrono::high_resolution_clock::now();
  e = y-x*betahat;
  arma::vec sigma = 1/sqrt(sum(x%x,1)/n);
  arma::vec psi = sigma%(as<arma::vec>(dnorm(e%sigma)));
  arma::mat xw = as<arma::mat>(sweep(x,1,sqrt(psi),"*"));
  D1 = xw.t()*xw/n;
  int counter = 0;
  for(int i = 0; i < n; i++){
    arma::mat xi = x.rows(counter,counter+rep(i)-1);
    arma::mat M = arma::zeros(rep(i), rep(i));
    if(type=="exchangeable"){
      M = arma::ones(rep(i), rep(i));
      M.diag(0) = arma::zeros(rep(i));
    }
    else{
      M.diag(1) = arma::ones(rep(i)-1);
      M.diag(-1) = arma::ones(rep(i)-1);
    }
    xM.cols(counter,counter+rep(i)-1)=xi.t()*M;
    counter = counter+rep(i);
  }
  arma::mat xMw = as<arma::mat>(sweep(xM,2,sqrt(psi),"*"));
  D2 = xMw*xw/n;
  arma::mat D = D2*inv(D1);
  auto finish_prep1 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_prep1 = finish_prep1 - start_prep1;
  time_prep1 = elapsed_prep1.count();
  time = time+time_prep1;
  
  //calculate g
  arma::mat g = arma::zeros(n, p);
  counter = 0;
  auto start_prep2 = std::chrono::high_resolution_clock::now();
  for(int i = 0; i < n; i++){
    arma::mat xi = x.rows(counter, counter+rep(i)-1);
    arma::vec lossi = lossweight(e.subvec(counter,counter+rep(i)-1), tau);
    g.row(i) = ((xM.cols(counter,counter+rep(i)-1)-D*xi.t())*lossi).t();
    counter = counter+rep(i);
  }
  auto finish_prep2 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_prep2 = finish_prep2-start_prep2;
  time_prep2 = elapsed_prep2.count();
  time = time+time_prep2;
  
  //update lambda
  arma::vec lambdaini = arma::zeros(p), lambda = arma::zeros(p);
  int iteration = 0;
  double distance = 1;
  while((iteration<maxstep)&&(distance>eps)){
    arma::vec S = arma::zeros(p);
    arma::mat F = arma::zeros(p, p);
    auto start_iter = std::chrono::high_resolution_clock::now();
    arma::vec tmp = 1/(g*lambdaini+arma::ones(n));
    S = g.t()*tmp;
    arma::mat gw = arma::zeros(n, p);
    gw = as<arma::mat>(sweep(g,1,tmp,"*"));
    F = gw.t()*gw;
    lambda = lambdaini+inv(F)*S;
    auto finish_iter = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_iter = finish_iter-start_iter;
    time_iter = elapsed_iter.count();
    time = time+time_iter;
    distance = accu(abs(lambda-lambdaini));
    lambdaini = lambda;
    iteration = iteration+1;
  }
  
  //calculate weight
  arma::vec weight = arma::zeros(n);
  auto start_wei = std::chrono::high_resolution_clock::now();
  weight = 1/(n*(arma::ones(n)+g*lambda));
  auto finish_wei = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_wei = finish_wei-start_wei;
  time_wei = elapsed_wei.count();
  time = time+time_wei;
  
  //return final value
  arma::vec final = arma::zeros(n+2);
  final.subvec(0,n-1) = weight;
  final(n) = time;
  final(n+1) = iteration;
  return final;
  
}

// Compute the weights
// [[Rcpp::export]]
Rcpp::List WS(arma::mat x, arma::vec y, arma::vec rep, double tau, arma::vec betahat, bool intercept, String corrtype = "ar1", double epsw = 1e-06, double maxstepw = 100){
  
  Rcpp::List final;
  int n = rep.size(), nsum = x.n_rows;
  arma::vec weight = arma::zeros(n+2);
  weight = WeightCPP(x, y, rep, tau, betahat, epsw, maxstepw, intercept, corrtype);
  arma::vec weightnew = arma::zeros(nsum);
  int counter = 0;
  for(int j=0;j<n;j++){
    weightnew.subvec(counter, (counter+rep(j)-1)) = weight(j)*arma::ones(rep(j));
    counter = counter+rep(j);
  }
  final = List::create(Named("Weight") = weightnew, Named("Time_W") = weight(n));
  return final;
  
}

// Compute the WQR estimator
// [[Rcpp::export]]
Rcpp::List WQRADMM(arma::mat x, arma::vec y, arma::vec rep, double tau, bool intercept, String esttype, bool warmstart = true, String corrtype = "ar1", double rhoCQR = 1, double rhoWQR = 1, double eps = 1e-04, double epsw = 1e-06, int maxstep = 5000, int maxstepw = 100){
  
  Environment base("package:base");
  Function sweep = base["sweep"];
  
  int p = x.n_cols, n = rep.size(), nsum = x.n_rows;
  if(intercept){
    p = p+1;
  }
  arma::vec result_CQR = arma::zeros(p+2);
  arma::vec result_WQR = arma::zeros(p+2);
  Rcpp::List final;
  if(esttype == "CQR"){
    result_CQR = WQRCPP(x, y, tau, rhoCQR, eps, maxstep, intercept, warmstart);
    final = List::create(Named("Estimation_CQR") = result_CQR.subvec(0,p-1),Named("Iteration_CQR") = result_CQR(p+1), Named("Time_CQR") = result_CQR(p));
  }
  if(esttype == "WQR"){
    result_CQR = WQRCPP(x, y, tau, rhoCQR, eps, maxstep, intercept, warmstart);
    arma::vec weight = arma::zeros(n+2);
    weight = WeightCPP(x, y, rep, tau, result_CQR.subvec(0, p-1), epsw, maxstepw, intercept, corrtype);
    arma::vec weightnew = arma::zeros(nsum);
    int counter = 0;
    for(int j=0;j<n;j++){
      weightnew.subvec(counter, (counter+rep(j)-1)) = weight(j)*arma::ones(rep(j));
      counter = counter+rep(j);
    }
    arma::mat xw = as<arma::mat>(sweep(x,1,weightnew*n,"*"));
    arma::vec yw = (weightnew*n)%y;
    xw.insert_cols(0, weightnew*n); 
    result_WQR = WQRCPP(xw, yw, tau, rhoWQR, eps, maxstep, FALSE, warmstart);
    double Tsum = 0;
    Tsum = result_CQR(p)+weight(n)+result_WQR(p);
    final = List::create(Named("Estimation_CQR") = result_CQR.subvec(0,p-1), Named("Iteration_CQR") = result_CQR(p+1), Named("Time_CQR") = result_CQR(p), Named("Estimation_WQR") = result_WQR.subvec(0,p-1), Named("Iteration_WQR") = result_WQR(p+1), Named("Time_WQR") = result_WQR(p), Named("Time_total") = Tsum);
  }
  return final;
}


//The parallel WQR-ADMM algorithm
arma::vec paraWQRCPP(arma::mat x, arma::vec y, arma::vec rep, int K, double tau, double rho, double eps, int maxstep, bool intercept, bool warmstart){
  
  int nsum = x.n_rows;
  if(intercept){
    x.insert_cols(0, arma::ones(nsum));
  }
  int p = x.n_cols, n = rep.size(), nk = n/K;
  arma::vec betaini = arma::zeros(p), vini = arma::zeros(nsum); 
  arma::vec beta = arma::zeros(p), r = arma::zeros(nsum), v = arma::zeros(nsum); 
  arma::mat betakini = arma::zeros(p,K), betak = arma::zeros(p,K), ukini = arma::zeros(p,K), uk = arma::zeros(p,K);
  arma::vec xbeta = arma::zeros(nsum);
  double time = 0, time_warm = 0, max_prep = 0, time_reduce = 0, max_map = 0; 
  arma::vec time_map = arma::zeros(K);
  arma::vec final = arma::zeros(p+2);
  
  int counter = 0;
  arma::vec num = arma::zeros(K);
  if(warmstart){
    arma::vec z = arma::zeros(p);
    auto start_warm = std::chrono::high_resolution_clock::now();
    z = inv(x.t()*x)*x.t()*y;
    auto finish_warm = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_warm = finish_warm - start_warm;
    time_warm = elapsed_warm.count();
    time = time + time_warm;
    for(int k = 0; k < K; k++){
      num(k) = accu(rep.subvec(nk*k,nk*k+nk-1));
      arma::mat xk = x.rows(counter,counter+num(k)-1);
      betakini.col(k) = z;
      xbeta.subvec(counter,counter+num(k)-1) = xk*betakini.col(k);
      counter = counter+num(k);
    }
  }

  counter = 0;
  arma::cube dat = arma::zeros<arma::cube>(p,p,K);
  for(int k = 0; k < K; k++){
    num(k) = accu(rep.subvec(nk*k,nk*k+nk-1));
    arma::mat tmp = arma::zeros(p,p), xk = x.rows(counter,counter+num(k)-1);
    auto start_prep = std::chrono::high_resolution_clock::now();
    if(num(k) > p) 
      tmp = inv(xk.t()*xk+arma::eye(p,p));
    else tmp = arma::eye(p,p)-xk.t()*inv(xk*xk.t()+arma::eye(num(k),num(k)))*xk;
    auto finish_prep = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_prep = finish_prep - start_prep;
    if(elapsed_prep.count() > max_prep) max_prep = elapsed_prep.count();
    dat.slice(k) = tmp;
    counter = counter+num(k);
  }
  time = time+max_prep;
  
  int iteration = 0;
  double distance = 1;
  
  while((distance>eps)&&(iteration<maxstep)){
    
    //update beta
    auto start_reduce = std::chrono::high_resolution_clock::now();
    beta = mean(betakini,1)+mean(ukini,1)/rho;
    auto finish_reduce = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_reduce = finish_reduce - start_reduce;
    //time cost for beta-update in each iteration
    time_reduce = elapsed_reduce.count();   
    time = time+time_reduce;
    
    counter = 0;
    //update r, betak, uk, v
    for(int k=0;k<K;k++){
      
      arma::vec yk = y.subvec(counter,counter+num(k)-1), vk = vini.subvec(counter,counter+num(k)-1);
      arma::mat xk = x.rows(counter,counter+num(k)-1);
      
      auto start_map = std::chrono::high_resolution_clock::now();
      //update r
      r.subvec(counter,counter+num(k)-1) = shrinkcpp(vk/rho+yk-xbeta.subvec(counter,counter+num(k)-1)-0.5*(2*tau-1)/rho, 0.5*arma::ones<arma::vec>(num(k))/rho);
      //update betak
      betak.col(k) = dat.slice(k)*(xk.t()*(yk-r.subvec(counter,counter+num(k)-1)+vk/rho)+beta-ukini.col(k)/rho);
      //update uk and v
      xbeta.subvec(counter,counter+num(k)-1) = xk*betak.col(k);
      uk.col(k) = ukini.col(k)+rho*(betak.col(k)-beta);
      v.subvec(counter,counter+num(k)-1) = vini.subvec(counter,counter+num(k)-1)+rho*(yk-xbeta.subvec(counter,counter+num(k)-1)-r.subvec(counter,counter+num(k)-1));
      auto finish_map = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> elapsed_map = finish_map - start_map;
      time_map(k) = elapsed_map.count();
      counter = counter+num(k);
    
    }
    //time cost for local computation in each iteration
    max_map = time_map(time_map.index_max());
    time = time+max_map;
    if(accu(abs(betaini))==0)
      distance = 1;
    else distance = accu(abs(beta-betaini))/accu(abs(betaini));
    betaini = beta, betakini = betak, ukini = uk, vini = v;
    iteration = iteration + 1;
    
  }
 
  final.subvec(0,p-1) = beta;
  final(p) = time;
  final(p+1) = iteration;
  return final;
  
}

// The adjusted weight
arma::vec paraWeightCPP(arma::mat x, arma::vec y, arma::vec rep, int K, double tau, arma::vec betahat, double eps, double maxstep, bool intercept, String type){
  
  Environment stats("package:stats");
  Function dnorm = stats["dnorm"];
  Environment base("package:base");
  Function sweep = base["sweep"];
  
  int nsum = x.n_rows;
  if(intercept){
    x.insert_cols(0, arma::ones(nsum));
  }
  int p = x.n_cols, n = rep.size(), nk = n/K;
  double max_prep1 = 0, max_prep2 = 0, time_reduce1 = 0, time_reduce2 = 0, max_wei = 0, time = 0;
  arma::vec time_map = arma::zeros(K);
  
  //calculate D1 and D2
  arma::vec e = arma::zeros(nsum);
  arma::cube tmp1 = arma::zeros<arma::cube>(p,p,K), tmp2 = arma::zeros<arma::cube>(p,p,K);
  arma::mat xM = arma::zeros(p,nsum), D1 = arma::zeros(p,p), D2 = arma::zeros(p,p);
  int counter = 0;
  arma::vec num = arma::zeros(K);
  for(int k = 0; k < K; k++){
    num(k) = accu(rep.subvec(nk*k,nk*k+nk-1));
    arma::vec yk = y.subvec(counter, counter+num(k)-1);
    arma::mat xk = x.rows(counter, counter+num(k)-1);
    auto start_prep1=std::chrono::high_resolution_clock::now();
    e.subvec(counter, counter+num(k)-1) = yk-xk*betahat;
    arma::vec sigma = 1/sqrt(sum(xk%xk,1)/n);
    arma::vec psi = sigma%(as<arma::vec>(dnorm(e.subvec(counter,counter+num(k)-1)%sigma)));
    arma::mat xkw = as<arma::mat>(sweep(xk,1,sqrt(psi),"*"));
    tmp1.slice(k) = xkw.t()*xkw/n;
    double counteri = counter;
    for(int i = nk*k; i<(nk*k+nk); i++){
      arma::mat xi = x.rows(counteri,counteri+rep(i)-1);
      arma::mat M = arma::zeros(rep(i), rep(i));
      if(type=="exchangeable"){
        M = arma::ones(rep(i), rep(i));
        M.diag(0) = arma::zeros(rep(i));
      }
      else{
        M.diag(1) = arma::ones(rep(i)-1);
        M.diag(-1) = arma::ones(rep(i)-1);
      }
      xM.cols(counteri,counteri+rep(i)-1) = xi.t()*M;
      counteri = counteri+rep(i);
    }
    arma::mat xMw = as<arma::mat>(sweep(xM.cols(counter, counter+num(k)-1),2,sqrt(psi),"*"));
    tmp2.slice(k) = xMw*xkw/n;
    auto finish_prep1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_prep1 = finish_prep1 - start_prep1;
    if(elapsed_prep1.count()>max_prep1) max_prep1 = elapsed_prep1.count();
    counter = counter+num(k);
  }
  time = time+max_prep1;
  
  //sum D1 and D2 
  auto start_reduce1 = std::chrono::high_resolution_clock::now();
  for(int k = 0; k < K; k++){
    D1 = D1+tmp1.slice(k);
    D2 = D2+tmp2.slice(k);
  }
  arma::mat D = D2*inv(D1);
  auto finish_reduce1 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_reduce1 = finish_reduce1-start_reduce1;
  time_reduce1 = elapsed_reduce1.count();
  time = time+time_reduce1;
  
  //calculate g
  arma::mat g = arma::zeros(n, p);
  counter = 0;
  for(int k = 0; k < K; k++){
    auto start_prep2 = std::chrono::high_resolution_clock::now();
    for(int i = nk*k; i < nk*k+nk; i++){
      arma::mat xi = x.rows(counter, counter+rep(i)-1);
      arma::vec lossi = lossweight(e.subvec(counter,counter+rep(i)-1), tau);
      g.row(i) = ((xM.cols(counter,counter+rep(i)-1)-D*xi.t())*lossi).t();
      counter=counter+rep(i);
    }
    auto finish_prep2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_prep2 = finish_prep2-start_prep2;
    if(elapsed_prep2.count()>max_prep2) max_prep2 = elapsed_prep2.count();
  }
  time = time+max_prep2;
  
  //update lambda
  arma::vec lambdaini = arma::zeros(p), lambda = arma::zeros(p);
  int iteration = 0;
  double distance = 1;
  
  while((iteration < maxstep)&&(distance > eps)){
    arma::vec S = arma::zeros(p);
    arma::mat Sk = arma::zeros(p, K);
    arma::mat F = arma::zeros(p, p);
    arma::cube Fk = arma::zeros<arma::cube>(p,p,K);
    for(int k = 0; k < K; k++){
      arma::mat gk = g.rows(nk*k, nk*k+nk-1);
      auto start_map = std::chrono::high_resolution_clock::now();
      arma::vec tmp3 = 1/(gk*lambdaini+arma::ones(nk));
      Sk.col(k) = gk.t()*tmp3;
      arma::mat gw = as<arma::mat>(sweep(gk,1,tmp3,"*"));
      Fk.slice(k) = gw.t()*gw;
      auto finish_map = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> elapsed_map = finish_map - start_map;
      time_map(k) = elapsed_map.count();
    }
    time = time+time_map(time_map.index_max());
    auto start_reduce2 = std::chrono::high_resolution_clock::now();
    for(int k = 0; k < K; k++){
      F = F+Fk.slice(k);
      S = S+Sk.col(k);
    }
    lambda = lambdaini+inv(F)*S;
    auto finish_reduce2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_reduce2 = finish_reduce2-start_reduce2;
    time_reduce2 = elapsed_reduce2.count();
    time = time+time_reduce2;
    distance = sum(abs(lambda-lambdaini));
    lambdaini = lambda;
    iteration = iteration+1;
  }
  
  //calculate weight
  arma::vec weight = arma::zeros(n);
  for(int k = 0; k < K; k++){
    arma::mat gk = g.rows(nk*k, nk*k+nk-1);
    auto start_wei = std::chrono::high_resolution_clock::now();
    weight.subvec(nk*k, nk*k+nk-1) = 1/(n*(arma::ones(nk)+gk*lambda));
    auto finish_wei = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_wei = finish_wei-start_wei;
    if(elapsed_wei.count()>max_wei) max_wei = elapsed_wei.count();
  }
  time = time+max_wei;
  
  arma::vec final = arma::zeros(n+2);
  final.subvec(0,n-1) = weight;
  final(n) = time;
  final(n+1) = iteration;
  return final;
  
}

// Compute the WQR estimator in parallel
// [[Rcpp::export]]
Rcpp::List paraWQRADMM(arma::mat x, arma::vec y, int K, arma::vec rep, double tau, bool intercept, String esttype, bool warmstart = true, String corrtype = "ar1", double rhoCQR = 1, double rhoWQR = 1, double eps = 1e-04, double epsw = 1e-06, int maxstep = 5000, int maxstepw = 100){
  
  Environment base("package:base");
  Function sweep = base["sweep"];
  
  int p = x.n_cols, n = rep.size(), nsum = x.n_rows;
  if(intercept){
    p = p+1;
  }
  arma::vec pararesult_CQR = arma::zeros(p+2);
  arma::vec pararesult_WQR = arma::zeros(p+2);
  Rcpp::List final;
  if(esttype == "CQR"){
    pararesult_CQR = paraWQRCPP(x, y, rep, K, tau, rhoCQR, eps, maxstep, intercept, warmstart);
    final = List::create(Named("Estimation_CQR") = pararesult_CQR.subvec(0,p-1),Named("Iteration_CQR") = pararesult_CQR(p+1), Named("Time_CQR") = pararesult_CQR(p));
  }
  if(esttype == "WQR"){
    pararesult_CQR = paraWQRCPP(x, y, rep, K, tau, rhoCQR, eps, maxstep, intercept, warmstart);
    arma::vec paraweight = arma::zeros(n+2);
    paraweight = paraWeightCPP(x, y, rep, K, tau, pararesult_CQR.subvec(0, p-1), epsw, maxstepw, intercept, corrtype);
    arma::vec paraweightnew = arma::zeros(nsum);
    int counter = 0;
    for(int j=0;j<n;j++){
      paraweightnew.subvec(counter, (counter+rep(j)-1)) = paraweight(j)*arma::ones(rep(j));
      counter = counter+rep(j);
    }
    arma::mat xw = as<arma::mat>(sweep(x,1,paraweightnew*n,"*"));
    arma::vec yw = (paraweightnew*n)%y;
    xw.insert_cols(0, paraweightnew*n);
    pararesult_WQR = paraWQRCPP(xw, yw, rep, K, tau, rhoWQR, eps, maxstep, FALSE, warmstart);
    double Tsum = 0;
    Tsum = pararesult_CQR(p)+paraweight(n)+pararesult_WQR(p);
    final = List::create(Named("Estimation_CQR") = pararesult_CQR.subvec(0,p-1),Named("Iteration_CQR") = pararesult_CQR(p+1), Named("Time_CQR") = pararesult_CQR(p), Named("Estimation_WQR") = pararesult_WQR.subvec(0,p-1), Named("Iteration_WQR") = pararesult_WQR(p+1), Named("Time_WQR") = pararesult_WQR(p), Named("Time_total") = Tsum);
  }
  return final;
}











