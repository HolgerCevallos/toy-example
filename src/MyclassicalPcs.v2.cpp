#include <RcppArmadillo.h>

using namespace Rcpp ;
using namespace arma ;

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


Rcpp::List findpcs0(arma::mat x, arma::mat B, arma::rowvec mu, 
                   int Npc, double tol)
  {
  int iter = 0;
  int p = B.n_rows;
//  int h = hsub.n_rows;
  int h = x.n_rows;
  int q = B.n_cols;
  double sc0,sc,delta;
  arma::mat xcent;
  //arma::mat A(n,q); // 
  arma::mat aprov(q,h); // Note that this is the transposed matrix for convenience
  arma::mat Bprov=B.t(); // NUEVO!
  arma::mat xcentprov; // NUEVO!
  arma::mat Bprod; // NUEVO!
  arma::mat Aprod; // NUEVO!
  arma::mat xnorm;
  arma::mat Q, R; // NUEVO!
  arma::mat Bortho(p,q); // NUEVO!
  
      xcent=x.each_row() - mu;
      xcentprov = xcent.t(); // NUEVO!
      aprov = Bprov*xcentprov; // NUEVO! trans(xcent*B);
      xnorm = xcent - aprov.t()*Bprov;  // NUEVO! - trans(B*aprov);
      colvec sqresid_dist = sum((xnorm % xnorm),1);
      iter = 0;
      sc0 = sum(sqresid_dist) / h;
      delta = 1; // delta = 1-sc/sc0;
                // sc0 = sc;
            
      while(++iter < Npc && delta>tol )
        {
        // xcent = x.each_row() - mu;
//      Rcpp::Rcout << "The matrix is " << xcent << std::endl;
        try
          {
          Bprod = Bprov*Bprov.t(); // NUEVO! trans(B) * B; // NUEVO!
          for (int kloop =0; kloop<h; kloop++)
            {
//          arma::vec aprov = solve( trans(B) * B, trans(B) * hsubcent.col(kloop));
            aprov.col(kloop) = solve( Bprod , Bprov*xcentprov.col(kloop)); // NUEVO! trans(xcent.row(kloop)*B)); // NUEVO!
            } 
//          Rcpp::Rcout << "The matrix A is " << aprov << std::endl;
          Aprod = aprov * aprov.t(); // NUEVO! aprov * trans(aprov); // NUEVO!
          for(int jloop =0; jloop<p; jloop++) 
            {
            // B.row(jloop) = trans(solve( Aprod , aprov*xcent.col(jloop))); // NUEVO!
            Bprov.col(jloop) = solve( Aprod , aprov*xcent.col(jloop)); // NUEVO!
//          Rcpp::Rcout << "The product is " << mean(trans(B.row(jloop)*aprov)) << std::endl;
            // mu(jloop) = mean( x.col(jloop) - trans(B.row(jloop)*aprov)); 
            }
//          Rcpp::Rcout << "The matrix is " << B  << std::endl;
//          Rcpp::Rcout << "The center in the loop is " << mu << std::endl;      
          // xcent = x.each_row() - mu;
          xnorm = xcent - aprov.t()*Bprov;// NUEVO! -  trans(B*aprov);
          sqresid_dist = sum((xnorm % xnorm),1);
          sc = sum(sqresid_dist) / h; // why not mean??
          delta = 1-sc/sc0;
          sc0 = sc;
//          Rcpp::Rcout << "The scale in the loop is " << sc << std::endl;      
          }
        catch(std::runtime_error) 
          {	
          arma::vec eigval;
          arma::mat eigvec;
          // mu = mean(x);
          eig_sym(eigval, eigvec, cov(x));
          Bprov = trans(eigvec.cols(0,q-1)); // NUEVO!
          iter=Npc;
          }
        }
      B = Bprov.t(); // NUEVO!
      qr_econ(Q,R,B); // NUEVO!
      Bortho = Q; // NUEVO!
return Rcpp::List::create(Rcpp::Named("Bmat") = B, Rcpp::Named("Borthomat") = Bortho, Rcpp::Named("mu") = mu);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//


