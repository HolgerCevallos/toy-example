#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp ;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]

Rcpp::List findpcs(arma::mat x, arma::mat B, arma::rowvec mu, int h, double s, 
                   int N1, int N2, int Npc, double tol, bool fixed)
  {
  // int iter = 0, it=0; // NUEVO!
  int iter;
  int it=0; // NUEVO!
  int p = B.n_rows;
//  int h = hsub.n_rows;
  int n = x.n_rows;
  int q = B.n_cols;
  double sc0,sc,s0,delta1,delta2;
  arma::uvec x_index1, x_index2; // NUEVO!
  arma::mat xcent;
  arma::mat xcentprov; // NUEVO!
  arma::mat hsub, hsubcent, hsubcentprov; // NUEVO!
  arma::mat A(n,q); // NUEVO!
  arma::mat Aprov(q,n); // NUEVO!
  arma::mat aprov(q,h); // Note that this is the transposed matrix for convenience
  arma::mat amatrix; // NUEVO!
  arma::mat Bprov=B.t(); // NUEVO!
  arma::mat amatrixBprov; // NUEVO!
  arma::mat Bprod; // NUEVO!
  arma::mat aprod; // NUEVO!
  arma::mat hsubnorm;
  arma::mat xnorm;
  arma::mat Q, R; // NUEVO!
  arma::mat Bortho; // NUEVO!
  
//  Rcpp::Rcout << "B is " << B << std::endl; 
//  Rcpp::Rcout << "The center is " << mu << std::endl; 
  xcent=x.each_row() - mu; //centrar para llegar a calcular la resta  
  xcentprov = xcent.t(); // NUEVO!
  if ((N1>0) || fixed)
    {
    Aprov = Bprov*xcentprov; // NUEVO! A = xcent*B;//calculo de Aq
//    Rcpp::Rcout << "The scores are " << A << std::endl;
    }
  else
    {
    Bprod = Bprov * B; // NUEVO! Bprod = trans(B) * B; // NUEVO!
    for (int kloop = 0; kloop < n; kloop++)
      {
      //          arma::vec aprov = solve( trans(B) * B, trans(B) * hsubcent.col(kloop));
      Aprov.col(kloop) = solve( Bprod , Bprov*xcentprov.col(kloop));// NUEVO! , trans(xcent.row(kloop)*B)));
      } 
    }
  xnorm = xcent -  Aprov.t()*Bprov; // NUEVO! -  A*trans(B);
  colvec sqresid_dist = sum((xnorm % xnorm),1); 
//  Rcpp::Rcout << "The distances are " << trans(sqresid_dist) << std::endl;
  x_index2 = sort_index(sqresid_dist); // NUEVO!
  // x_index = sort_index(sqresid_dist); // NUEVO!
  delta1 = 1;
  s0 = s;
  
  while( ++it <= (N1+N2) && delta1>tol)
    {   
//  Make h subset
    //    ws<-rep(1,n)
    x_index1 = x_index2; // NUEVO!
    // x_index1 = sort_index(sqresid_dist); // NUEVO! // NUEVO!
    hsub = x.rows(x_index1.head(h)); // NUEVO!
    mu = mean(hsub); //creo que aquí calcula m, en función pesos. 
//    Rcpp::Rcout << "The subset is " << trans(sort(x_index.head(h)+1)) << std::endl;
//    Rcpp::Rcout << "The center is " << mu << std::endl;
//    Rcpp::Rcout << "The condition is " <<(!fixed & (it>N1)) << std::endl;
    if(!fixed & (it>N1))    
      {
      iter = 0;
//      sc0 = 1e10;
//      sc0 = std::numeric_limits<double>::infinity();
      if (s == 0) {sc0 = sqrt( sum(sqresid_dist(x_index1.head(h))) / h );}
      else {sc0 = s;}
      delta2 =1; // delta = 1-sc/sc0;
                // sc0 = sc;
      hsubcent = hsub.each_row() - mu; // NUEVO!       
      while(++iter < Npc && delta2>tol )
        {
        // hsubcent = hsub.each_row() - mu; // NUEVO!
        hsubcentprov = hsubcent.t(); // NUEVO! 
//      Rcpp::Rcout << "The matrix is " << hsubcent << std::endl;
        try
          {
          Bprod = Bprov * Bprov.t(); // NUEVO! trans(B) * B; // NUEVO!
          for (int kloop =0; kloop<h; kloop++)
            {
//          arma::vec aprov = solve( trans(B) * B, trans(B) * hsubcent.col(kloop));
            aprov.col(kloop) = solve( Bprod , Bprov*hsubcentprov.col(kloop)); // NUEVO! trans(hsubcent.row(kloop)*B)); // NUEVO!
            } 
//          Rcpp::Rcout << "The matrix A is " << aprov << std::endl;
          amatrix = aprov.t(); // NUEVO!
          aprod = aprov*amatrix; // NUEVO! aprov * trans(aprov); // NUEVO!
          for(int jloop =0; jloop<p; jloop++) 
            {
            // B.row(jloop) = trans(solve( Aprod , aprov*hsubcent.col(jloop))); // NUEVO!
            Bprov.col(jloop) = solve( aprod , aprov*hsubcent.col(jloop)); // NUEVO!
            //  Rcpp::Rcout << "The product is " << mean(trans(B.row(jloop)*aprov)) << std::endl;
            // mu(jloop) = mean( hsub.col(jloop) - amatrix*Bprov.col(jloop)); // NUEVO ! - trans(B.row(jloop)*aprov)); 
            }
//          Rcpp::Rcout << "The matrix is " << B  << std::endl;
//          Rcpp::Rcout << "The center in the loop is " << mu << std::endl;      
          amatrixBprov = amatrix*Bprov; // NUEVO!
          mu = mean(hsub - amatrixBprov); // NUEVO!
          hsubcent = hsub.each_row() - mu;
          hsubnorm = hsubcent - amatrixBprov; // NUEVO! -  trans(B*aprov);
          colvec sqresid_dist0 = sum((hsubnorm % hsubnorm),1);
          sc = sqrt( sum(sqresid_dist0) / h ); // why not mean??
          delta2 = 1-sc/sc0;
          sc0 = sc;
//          Rcpp::Rcout << "The scale in the loop is " << sc << std::endl;      
          }
        catch(std::runtime_error) 
          {	
          arma::vec eigval;
          arma::mat eigvec;
          mu = mean(hsub);
          eig_sym(eigval, eigvec, cov(hsub));
          Bprov = trans(eigvec.cols(0,q-1)); // NUEVO!
          iter=Npc;
          }
        }
      xcent=x.each_row() - mu;
      xcentprov = xcent.t(); // NUEVO!
      Bprod = Bprov * Bprov.t();// NUEVO! trans(B) * B; // NUEVO!
      for (int kloop = 0; kloop < n; kloop++)
        {
        //          arma::vec aprov = solve( trans(B) * B, trans(B) * hsubcent.col(kloop));
        Aprov.col(kloop) = solve( Bprod , Bprov*xcentprov.col(kloop)); // NUEVO!
        // A.row(kloop) = trans(solve( Bprod , trans(xcent.row(kloop)*B))); // NUEVO!
        } 
      }
    else
      {
      xcent=x.each_row() - mu;
      xcentprov = xcent.t(); // NUEVO!
      Aprov = Bprov*xcentprov; // NUEVO! A = xcent*B;
//      Rcpp::Rcout << "xcent is " << xcent.row(0)  << std::endl;
//      Rcpp::Rcout << "xcent is " << xcent  << std::endl;
//      Rcpp::Rcout << "B is " << B  << std::endl;
//      Rcpp::Rcout << "A is " << A.row(1)  << std::endl;
     // A = (x.each_row() - mu) * B;
//      Rcpp::Rcout << "The first row is " << A  << std::endl;
      }
    xnorm = xcent -  Aprov.t()*Bprov; // NUEVO! -  A * trans(B);
//    Rcpp::Rcout << "xnorm is " << xnorm.row(0)  << std::endl;
    sqresid_dist = sum((xnorm % xnorm),1);
//    Rcpp::Rcout << "dist is " << sqresid_dist  << std::endl;
    // x_index = sort_index(sqresid_dist); // NUEVO!
    x_index2 = sort_index(sqresid_dist); // NUEVO!
    if (N1==0 && s!=0)
      {
//      Rcpp::Rcout << "B is " << B << std::endl;   
      // x_index2 = sort_index(sqresid_dist); // NUEVO!
      s = sqrt( sum(sqresid_dist(x_index2.head(h))) / h ); // NUEVO!
      delta1 = 1-s/s0;
      s0 = s;
//      Rcpp::Rcout << "The scale is " << s << std::endl;      
      }
    }
  if(s==0)
    {
    // x_index = sort_index(sqresid_dist); // NUEVO!
    s = sqrt( sum(sqresid_dist(x_index2.head(h))) / h ); // NUEVO! 
    }
  x_index1 = sort(x_index1.head(h)+1); // NUEVO!
  // x_index2 = x_index2.head(h)+1; // NUEVO!
//  return(B);

if(fixed){ // NUEVO!
  Bortho = B; // NUEVO!
} 
else // NUEVO!
{ 
  B = Bprov.t(); // NUEVO!
  qr_econ(Q,R,B); // NUEVO!
  Bortho = Q; // NUEVO!
}

A = (x.each_row() - mu)*Bortho; // NUEVO!

return Rcpp::List::create(Rcpp::Named("Bmat") = B, Rcpp::Named("Borthomat") = Bortho, Rcpp::Named("mu") = mu, Rcpp::Named("scale") = s,
                          Rcpp::Named("index1") = x_index1, // NUEVO!
                          Rcpp::Named("Amat") = A); // NUEVO!
                            //, Rcpp::Named("index2") = x_index2); // NUEVO!
// return Rcpp::List::create(Rcpp::Named("Bmat") = mu);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//


