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


// [[Rcpp::export]]
Rcpp::List get_sequence(Rcpp::StringVector chrom,
                  Rcpp::IntegerVector start,
                  Rcpp::IntegerVector end,
                  Rcpp::StringVector strand,
                  Rcpp::List genome)
{
  int interval_count = chrom.length();
  //std::cout << "Number of intervals: " << interval_count << "\n";
  //int chromosome_count = genome.length();
  //std::cout << "Number of chromosomes: " << chromosome_count << "\n";
  //Rcpp::StringVector chromosome_names = genome.names();
  //std::cout << "Chromosome names: " << chromosome_names[2] << "\n";
  //std::string ChrID = "MtDNA";
  //std::string chromosome = genome[ChrID];
  //std::cout << "Chromosome 1 length: " << ChrID << "    "<<chromosome << "\n";
  //std::cout << "mtDNA: " << Rcpp::String(genome["MtDNA"])<< "\n";

  //for(int i = 0; i<= length; i++){
  //  Rcpp::List::create(Rcpp::Named("DNA") = genome(i).substr(start,end-start+1));
  //}
//
//   Rcpp::CharacterVector results1 = Rcpp::CharacterVector::create(Rcpp::Named("fu3") = chromosome.substr(7,1),
//                                                                Rcpp::Named("fu2") = chromosome.substr(8,1),
//                                                                Rcpp::Named("fu1") = chromosome.substr(9,1));
//   Rcpp::CharacterVector results2 = Rcpp::CharacterVector::create(Rcpp::Named("fu3") = chromosome.substr(17,1),
//                                                                Rcpp::Named("fu2") = chromosome.substr(18,1),
//                                                                Rcpp::Named("fu1") = chromosome.substr(19,1));
  // Rcpp::List new_list(4);
  // new_list.push_back(results1);
  // new_list.push_back(results2);
  // new_list.push_back(results1);
  // new_list.push_back(results2);

  Rcpp::List new_list(interval_count);
  // Rcpp::StringVector Results_Interval = Rcpp::StringVector(10);
  // Rcpp::StringVector Results_With_flanks = Rcpp::StringVector(10);
  // Rcpp::StringVector Results_Five = Rcpp::StringVector(10);
  // Rcpp::StringVector Results_Three = Rcpp::StringVector(10);


  //, "with_flanks","five", "three");
  std::string full_chromosome;
  std::string new_chromosome = "";
  std::string old_chromosome = "NA";
  int new_start;
  int length;
  std::string new_strand = "";
  std::string new_interval = "";
  for(int i = 0; i<interval_count; i++){
    new_chromosome = chrom[i];
    if(new_chromosome != old_chromosome){
      std::cout << "New" << new_chromosome << old_chromosome << "\n";
      std::cout << full_chromosome.length() << " " << start[i] << " " << end[i] << "\n";
      full_chromosome = Rcpp::as< std::string >(genome[new_chromosome]);
      std::cout << full_chromosome.length() << " " << chrom[i] << "\n";
      old_chromosome = new_chromosome;
    }
    //chromosome_name = chrom[i];
    //std::cout << "Chromosome names: " << new_chromosome << "\n";
    //std::string full_chromosome = genome[new_chromosome];
    //std::cout << "Substring: " << full_chromosome.substr(19, 3) << "\n";
    //std::cout << full_chromosome.length() << " " << start[i] << " " << end[i] << "\n";
    new_start = start[i]-1-10;
    if((new_start<0) || (new_start+length>full_chromosome.length())){
      new_interval = "NNNNNNNNNNNNNNNNNNNN"
    }
    length = end[i] - new_start + 10;
    new_strand = strand[i];
    new_interval = full_chromosome.substr(new_start, length);
    if(new_strand == "-"){
      new_interval = reverse_complement(new_interval);
    }

//
//     new_list[i+1] =  Rcpp::StringVector::create(new_interval.substr(10, new_interval.length()-20),
//                                                 new_interval,
//                                                 new_interval.substr(10,1),
//                                                 new_interval.substr(new_interval.length()-11,1));
//
//
//
//     Results_Interval[i] = Rcpp::StringVector::create(new_interval.substr(10, new_interval.length()-20));
//     Results_With_flanks[i] = new_interval;
//     Results_Five[i] = new_interval.substr(10,1);
//     Results_Three[i] = new_interval.substr(new_interval.length()-11,1);


    new_list[i] = Rcpp::StringVector::create(Rcpp::Named("interval") = new_interval.substr(10, new_interval.length()-20),
                                              Rcpp::Named("with_flanks") = new_interval,
                                              Rcpp::Named("fu3") = new_interval.substr(7,1),
                                              Rcpp::Named("fu2") = new_interval.substr(8,1),
                                              Rcpp::Named("fu1") = new_interval.substr(9,1),
                                              Rcpp::Named("five") = new_interval.substr(10,1),
                                              Rcpp::Named("fd1") = new_interval.substr(11,1),
                                              Rcpp::Named("fd2") = new_interval.substr(12,1),
                                              Rcpp::Named("fd3") = new_interval.substr(13,1),
                                              Rcpp::Named("tu3") = new_interval.substr(new_interval.length()-14,1),
                                              Rcpp::Named("tu2") = new_interval.substr(new_interval.length()-13,1),
                                              Rcpp::Named("tu1") = new_interval.substr(new_interval.length()-12,1),
                                              Rcpp::Named("three") = new_interval.substr(new_interval.length()-11,1),
                                              Rcpp::Named("td1") = new_interval.substr(new_interval.length()-10,1),
                                              Rcpp::Named("td2") = new_interval.substr(new_interval.length()-9,1),
                                              Rcpp::Named("td3") = new_interval.substr(new_interval.length()-8,1));
    //std::cout << "Chromosome 1 length: " << genome[z] << "\n";
    //full_sequence = &::String("test"); //&genome["I"];
    //std::cout << "Chromosome 1" << full_sequence << "\n";
    //std::cout << "Chromosome 1 subset: " << std::string(); //.substr(start[i], end[i]-start[i]+1); << "\n";
    //full_sequence = genome[chrom[i]]; //.substr(start[i], end[i]-start[i]+1);
    //new_list[i] = new_interval;
  }
  return new_list;
  // Rcpp::Named("five") = chromosome.substr(10,1),
  // Rcpp::Named("fd1") = chromosome.substr(11,1),
  // Rcpp::Named("fd2") = chromosome.substr(12,1),
  // Rcpp::Named("fd3") = chromosome.substr(13,1),
  // Rcpp::Named("tu3") = chromosome.substr(chromosome.length()-14,1),
  // Rcpp::Named("tu2") = chromosome.substr(chromosome.length()-13,1),
  // Rcpp::Named("tu1") = chromosome.substr(chromosome.length()-12,1),
  // Rcpp::Named("three") = chromosome.substr(chromosome.length()-11,1),
  // Rcpp::Named("td1") = chromosome.substr(chromosome.length()-10,1),
  // Rcpp::Named("td2") = chromosome.substr(chromosome.length()-9,1),
  //Rcpp::Named("td3") = chromosome.substr(chromosome.length()-8,1));
  // return Rcpp::List::create(Rcpp::Named("Interval") = Results_Interval,
  //                           Rcpp::Named("with_flanks") = Results_With_flanks,
  //                           Rcpp::Named("five") = Results_Five,
  //                           Rcpp::Named("three") = Results_Three);

  //return chromosome.substr(5,10);
  //return Rcpp::List::create(interval_count, chromosome_count, chromosome_names);
  //std::string sequence = DNAseq.substr(start-11,end-start+1);
  //  if(strand=="-"){
  //    sequence = reverse_complement(sequence);

  //  return Rcpp::List::create(Rcpp::Named("DNAseq") = sequence);
  //std::cout << DNAseq.substr(4,6) << "\n";
  /*
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
  */
}




// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

///*** R
//timesTwo(42)
//*/
