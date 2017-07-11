#include <Rcpp.h>
#include <string>     // std::string, std::stoi
#include <iostream>
#include <algorithm>
#include <sstream>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List CalculateOffset(int primaryPosition,
                           Rcpp::IntegerVector secondaryPosition,
                           Rcpp::IntegerVector secondaryWidth,
                           Rcpp::StringVector secondaryChromosome,
                           const int& maxOffset){
  // int secondaryLength = secondaryPosition.length();
  // Rcpp::IntegerVector widths = Rcpp::IntegerVector::create();
  // std::cout << secondaryLength << "\n";
  Rcpp::IntegerVector offsets = secondaryPosition - primaryPosition;
  Rcpp::LogicalVector offsetMatch = Rcpp::abs(offsets) <= maxOffset;
  //
  //   offset = secondaryPosition[i] - primaryPosition;
  //   if(Rcpp::as<std::string>(secondaryChromosome[i]) == primaryChromosome && std::abs(offset) <= maxOffset){
  //     offsets.push_back(offset);
  //     widths.push_back(secondaryWidth[i]);
  //     chromosomes.push_back(secondaryChromosome[i]);

  // std::cout << "Offsets" << offsets << "offsetMatch" << offsetMatch << "\n";
  return Rcpp::List::create(Named("offsets") = offsets[offsetMatch],
                            Named("widths") = secondaryWidth[offsetMatch],
                            Named("chromosomes") = secondaryChromosome[offsetMatch]);
}




// [[Rcpp::export]]
Rcpp::NumericVector calculate_offset(const NumericMatrix& x, const NumericVector& A, const NumericVector& B) {
  int x_row = x.nrow();
  Rcpp::NumericVector results = Rcpp::NumericVector(x_row);
  // int res;
  for(int i = 0; i < x_row; i++){
    results(i) = B(x(i,1)-1) - A(x(i,0)-1);
    // std::cout << i << " " << x(i,0) << " " << x(i,1) << " " << A(x(i,0)-1) << " " << B(x(i,1)-1) << "\n";
  }
  // std::cout << x_row << " "  << " " << results.length() << "\n";
  // std::cout << A(x(0,0)) << " " << B(x(0,1)) << "\n";
  // results(1) = x(1,4);
return results;
}







// [[Rcpp::export]]
int min_index(NumericVector x) {
  // Rcpp supports STL-style iterators
  NumericVector::iterator it = std::min_element(x.begin(), x.end());
  // we want the value so dereference
  return it - x.begin();
}


// [[Rcpp::export]]
Rcpp::NumericVector calculate_distance2(Rcpp::NumericVector A, Rcpp::NumericVector B){
  int length = A.length();
  Rcpp::NumericVector results = Rcpp::NumericVector(length);
  for(int i = 0; i < length; i++){
    //std::cout << i << " "<< A(i) << '\n';
    // results(i) = A(min_index(B-A(i)));
    Rcpp::NumericVector first = B-A(i);
    int second = min_index(first);
    //std::cout << "first" << first << "second" << second << '\n';
    results(i) = A(second);
  }
  //int results = A(min(A)/13);
  return results;
}





// [[Rcpp::export]]
Rcpp::NumericVector calculate_distance(Rcpp::NumericVector A, Rcpp::NumericVector B){
//int calculate_distance(Rcpp::NumericVector A, Rcpp::NumericVector B){
  //int results = min(A);
  Rcpp::NumericVector results = 5-A;
  return results;
  //return A-B;
//   int A_length = A.length();
//   int B_length = B.length();
//   Rcpp::NumericVector results = Rcpp::NumericVector(A_length);
//   int num;
//   for(int i = 0; i < A_length; i++){
//     int num = A(i);
//     int l = 0;
//     int h = B_length;
// //    while(num)

  // }

  // return results;
}


//
// // [[Rcpp::export]]
// Rcpp::LogicalVector filter_MD_tags3(std::string strand, Rcpp::StringVector MD){
//   //std::cout << strand << NM << MD << "\n";
//   int MD_length = MD.length();
//   //std::string value;
//   int MD_value;
//   std::string MD_string;
//   if (strand == "+"){
//     MD_string = Rcpp::as< std::string >(MD(0));
//     //std::cout << "+MD_String" << MD_string << "\n";
//     std::stringstream(MD_string) >> MD_value;
//     //std::cout << "+MD_Value" << MD_value << "\n";
//   }
//   if (strand == "-"){
//     MD_string = Rcpp::as< std::string >(MD(MD_length-1));
//     //std::cout << "-MD_String" << MD_string << "\n";
//     std::stringstream(MD_string) >> MD_value;
//     //std::cout << "-MD_Value" << MD_value << "\n";
//   }
//   if (MD_value >= 22){
//     return Rcpp::LogicalVector(true);
//   }
//   return Rcpp::LogicalVector(false);
// }
//
//
//
//
//
//
// // [[Rcpp::export]]
// Rcpp::LogicalVector filter_MD_tags(std::string strand, int NM, Rcpp::StringVector MD){
//   if (NM == 0){
//     return Rcpp::LogicalVector(true);
//   }
//   //std::cout << strand << NM << MD << "\n";
//   int MD_length = MD.length();
//   //std::string value;
//   int MD_value;
//   std::string MD_string;
//   if (strand == "+"){
//     MD_string = Rcpp::as< std::string >(MD(0));
//     //std::cout << "+MD_String" << MD_string << "\n";
//     std::stringstream(MD_string) >> MD_value;
//     //std::cout << "+MD_Value" << MD_value << "\n";
//   }
//   if (strand == "-"){
//     MD_string = Rcpp::as< std::string >(MD(MD_length-1));
//     //std::cout << "-MD_String" << MD_string << "\n";
//     std::stringstream(MD_string) >> MD_value;
//     //std::cout << "-MD_Value" << MD_value << "\n";
//   }
//   if (MD_value >= 22){
//     return Rcpp::LogicalVector(true);
//   }
//   return Rcpp::LogicalVector(false);
// }

// // [[Rcpp::export]]
// Rcpp::LogicalVector filter_MD_tag(std::string strand, int NM, std::string MD){
//   if (NM == 0){
//     return Rcpp::LogicalVector(true);
//   }
//   //std::cout << strand << NM << MD << "\n";
//   int MD_length = MD.length();
//   std::string value;
//   int MD_value;
//   std::string MD_string;
//   if (strand == "+"){
//     int i = 0;
//     while (i < MD_length){
//       value = MD.substr(i,1);
//       if(value == "A" || value == "C" || value == "T" || value == "G" || value == "N"){
//         break;
//       }
//       i++;
//     }
//     MD_string = MD.substr(0,i);
//   }
//
//   if (strand == "-"){
//     int i = MD_length;
//     while (i > 0){
//       i--;
//       value = MD.substr(i,1);
//       if(value == "A" || value == "C" || value == "T" || value == "G" || value == "N"){
//         i = i+1;
//         break;
//       }
//       //std::cout << strand << MD << i << MD_length << "\n";
//       //printf("%d", i);
//     }
//     MD_string = MD.substr(i,MD_length-i);
//   }
//   std::stringstream(MD_string) >> MD_value;
//   if (MD_value >= 22){
//     return Rcpp::LogicalVector(true);
//   }
//   return Rcpp::LogicalVector(false);
// }

// [[Rcpp::export]]
std::string reverse_complement(std::string DNAseq)
{
  // Reverse the string
  reverse(DNAseq.begin(), DNAseq.end());

  // Loop over string and replace each base with complement
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
  // Return new sequence
  return DNAseq;
}


// [[Rcpp::export]]
List sequence_slice(std::string DNAseq)
{
  // Return the slices of the sequence in a named list
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
}

// [[Rcpp::export]]
Rcpp::List get_sequence(Rcpp::StringVector chrom,
                         Rcpp::IntegerVector start,
                         Rcpp::IntegerVector end,
                         Rcpp::StringVector strand,
                         Rcpp::List genome)
{
  // Get the number of intervals to be processed
  int interval_count = chrom.length();

  // Create the StringVectors to store the results
  Rcpp::StringVector vect_interval = Rcpp::StringVector(interval_count);
  Rcpp::StringVector vect_with_flanks = Rcpp::StringVector(interval_count);
  Rcpp::StringVector vect_fu3 = Rcpp::StringVector(interval_count);
  Rcpp::StringVector vect_fu2 = Rcpp::StringVector(interval_count);
  Rcpp::StringVector vect_fu1 = Rcpp::StringVector(interval_count);
  Rcpp::StringVector vect_five = Rcpp::StringVector(interval_count);
  Rcpp::StringVector vect_fd1 = Rcpp::StringVector(interval_count);
  Rcpp::StringVector vect_fd2 = Rcpp::StringVector(interval_count);
  Rcpp::StringVector vect_fd3 = Rcpp::StringVector(interval_count);
  Rcpp::StringVector vect_tu3 = Rcpp::StringVector(interval_count);
  Rcpp::StringVector vect_tu2 = Rcpp::StringVector(interval_count);
  Rcpp::StringVector vect_tu1 = Rcpp::StringVector(interval_count);
  Rcpp::StringVector vect_three = Rcpp::StringVector(interval_count);
  Rcpp::StringVector vect_td1 = Rcpp::StringVector(interval_count);
  Rcpp::StringVector vect_td2 = Rcpp::StringVector(interval_count);
  Rcpp::StringVector vect_td3 = Rcpp::StringVector(interval_count);

  // Variable for full chromosome sequence
  std::string full_chromosome;

  // Name of chromosome being processed
  std::string new_chromosome = "";

  // Name of last chromosome processed
  std::string old_chromosome = "NA";

  // Variables for start position, length, strand and sequence
  int new_start;
  int length;
  std::string new_strand = "";
  std::string new_interval = "";

  // Loop over each interval
  for(int i = 0; i<interval_count; i++){

    // Get the ID of the chromosome for this interval
    new_chromosome = chrom[i];

    // If the chromosome has changed, get the sequence for the new chromosome
    if(new_chromosome != old_chromosome){
      full_chromosome = Rcpp::as< std::string >(genome[new_chromosome]);
      old_chromosome = new_chromosome;
    }

    // Calculate new start position
    new_start = start[i]-1-10;

    // If the sequence goes off the end of the chromosome, return "N"
    if((new_start<0) || (new_start+length>full_chromosome.length())){
      new_interval = "NNNNNNNNNNNNNNNNNNNN";
    }

    // Calculate the length of the interval
    length = end[i] - new_start + 10;

    // Get the strand of the interval
    new_strand = strand[i];

    // Get the sequence of the interval
    new_interval = full_chromosome.substr(new_start, length);

    // If minus strand, reverse complement the sequence
    if(new_strand == "-"){
      new_interval = reverse_complement(new_interval);
    }

    // Store the results in the StringVectors at position "i"
    vect_interval[i] = new_interval.substr(10, new_interval.length()-20);
    vect_with_flanks[i] = new_interval;
    vect_fu3[i] = new_interval.substr(7,1);
    vect_fu2[i] = new_interval.substr(8,1);
    vect_fu1[i] = new_interval.substr(9,1);
    vect_five[i] = new_interval.substr(10,1);
    vect_fd1[i] = new_interval.substr(11,1);
    vect_fd2[i] = new_interval.substr(12,1);
    vect_fd3[i] = new_interval.substr(13,1);
    vect_tu3[i] = new_interval.substr(new_interval.length()-14,1);
    vect_tu2[i] = new_interval.substr(new_interval.length()-13,1);
    vect_tu1[i] = new_interval.substr(new_interval.length()-12,1);
    vect_three[i] = new_interval.substr(new_interval.length()-11,1);
    vect_td1[i] = new_interval.substr(new_interval.length()-10,1);
    vect_td2[i] = new_interval.substr(new_interval.length()-9,1);
    vect_td3[i] = new_interval.substr(new_interval.length()-8,1);
  }
  // Return a named list of the StringVectors
  return Rcpp::List::create(Rcpp::Named("interval") = vect_interval,
                            Rcpp::Named("with_flanks") = vect_with_flanks,
                            Rcpp::Named("fu3") = vect_fu3,
                            Rcpp::Named("fu2") = vect_fu2,
                            Rcpp::Named("fu1") = vect_fu1,
                            Rcpp::Named("five") = vect_five,
                            Rcpp::Named("fd1") = vect_fd1,
                            Rcpp::Named("fd2") = vect_fd2,
                            Rcpp::Named("fd3") = vect_fd3,
                            Rcpp::Named("tu3") = vect_tu3,
                            Rcpp::Named("tu2") = vect_tu2,
                            Rcpp::Named("tu1") = vect_tu1,
                            Rcpp::Named("three") = vect_three,
                            Rcpp::Named("td1") = vect_td1,
                            Rcpp::Named("td2") = vect_td2,
                            Rcpp::Named("td3") = vect_td3);
}
