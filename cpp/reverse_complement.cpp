#include <Rcpp.h>
#include <string>
#include <iostream>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
std::string reverse_complement(std::string DNAseq)
{
  //std::string DNAseq = "TGAGACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGC";
  reverse(DNAseq.begin(), DNAseq.end());
  for (std::size_t i = 0; i < DNAseq.length(); ++i){
    switch (DNAseq[i]){
    case 'A':
      DNAseq[i] = 'T';
      break;
    case 'C':
      DNAseq[i] = 'G';
      break;
    case 'G':
      DNAseq[i] = 'C';
      break;
    case 'T':
      DNAseq[i] = 'A';
      break;
    }
  }
  //std::cout << DNAseq; //<< std::endl;
  return DNAseq;
}


// [[Rcpp::export]]
List sequence_slice(std::string DNAseq)
{
  //std::cout << DNAseq.substr(4,6) << "\n";
  return Rcpp::List::create(Rcpp::Named("fu3") = DNAseq.substr(7,1),
                           Rcpp::Named("fu2") = DNAseq.substr(8,1),
                           Rcpp::Named("fu1") = DNAseq.substr(9,1),
                           Rcpp::Named("five") = DNAseq.substr(10,1),
                           Rcpp::Named("fd1") = DNAseq.substr(11,1),
                           Rcpp::Named("fd2") = DNAseq.substr(12,1),
                           Rcpp::Named("fd3") = DNAseq.substr(13,1),
                           Rcpp::Named("tu3") = DNAseq.substr(DNAseq.length()-14,1),
                           Rcpp::Named("tu2") = DNAseq.substr(DNAseq.length()-13,1),
                           Rcpp::Named("tu1") = DNAseq.substr(DNAseq.length()-12,1),
                           Rcpp::Named("three") = DNAseq.substr(DNAseq.length()-11,1),
                           Rcpp::Named("td1") = DNAseq.substr(DNAseq.length()-10,1),
                           Rcpp::Named("td2") = DNAseq.substr(DNAseq.length()-9,1),
                           Rcpp::Named("td3") = DNAseq.substr(DNAseq.length()-8,1));
  //std::cout << DNAseq; //<< std::endl;
  //return DNAseq;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

///*** R
//timesTwo(42)
//*/
