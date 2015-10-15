#include <string>
#include <iostream>
#include <fstream>

#include <cstdlib>
#include <cmath>
#include <cstdio>

#include "config.h"

double qualsMap[256];

// WARNING: this is the maximum size of a line for the input fastq file. Change
// it if larger lines need to be supported (e.g., extra long reads like PacBio).
#define MAX_BUFF_SIZE 4096
char buffer[MAX_BUFF_SIZE];

struct Read {
  Read() : header(""), bases(""), extra(""), qualities("") {}
  std::string header;
  std::string bases;
  std::string extra;
  std::string qualities;
};


Read nextRead(std::ifstream& is) {
  Read r;
  // header
  is.getline(buffer, MAX_BUFF_SIZE);
  if (!is.eof()) {
    r.header = buffer;

    // bases
    is.getline(buffer, MAX_BUFF_SIZE);
    r.bases = buffer;
  
    // extra
    is.getline(buffer, MAX_BUFF_SIZE);
    r.extra = buffer;

    // qualities
    is.getline(buffer, MAX_BUFF_SIZE);
    r.qualities = buffer;
  }
  return r;
}

void writeRead(const Read& r, std::ofstream& os) {
  os << r.header << '\n' << r.bases << '\n' << r.extra << '\n' << r.qualities << '\n';
}

void init() {
  for(size_t i = 0; i < 256; ++i) {
    qualsMap[i] = pow(10, -1.0 * ( (float)i / 10.0 ) );
  }
}

// transforms an error probability p into an encoded quality score Q using the
//    Q = 33 + floor(-10log10(p))
// equation for inverted phred scaled encoded using Sanger fastq format.
char probToScoreChar(double p) {
  // printf("%.15f %.15f\n",p, -10.0 * log10((double)(1.0-p)));
  char c = ('!' + round(-10.0 * log10((double)(1.0-p))));
  return c;
}

// Here is the 'core' function of the software that takes a qualities string
// and re-encode it to represent k-mer quality (BUG-FREE INEFFICIENT VERSION)
std::string reEncodeQuality(const std::string& q, size_t k) {
  size_t m = q.length();

  std::string encQ("");
  encQ += q[0];
  for (size_t i = 1; i < m-1; ++i) {
    int b = ( ((int)i-((int)k-1)) > 0) ? i-k+1 : 0;
    double p = 1.0;
    for (size_t j = b; j <= i; ++j) {
      //  std::cout << "[" << q[j] <<"] " << 1.0 - qualsMap[q[j] - '!'] << " * " << p ;
      p *= 1.0 - qualsMap[q[j] - '!'];
      //std::cout << " (" <<probToScoreChar(p) <<")"<< std::endl;
    }
    encQ += probToScoreChar(p);
  }
  encQ += q[m-1];
  return encQ;
}


// This procedure is more efficient, but does not work properly when there are
// 0 quality bases because the algorithm is not able to 'recover' from a 'fall'
// to a 0 probability (i.e., it then divides by zero getting -inf probability
// of correct base not matter what comes next)


/*std::string reEncodeQuality(const std::string& q, size_t k) {
  size_t m = q.size();
  std::string encQ(m,'~');

  // the first quality is unchanged
  encQ[0] = q[0];
  double currP = 1.0 - qualsMap[q[0] - 33];
  for (size_t i = 1; i < k; ++i) {
    std::cout << currP << " " << 1.0 - qualsMap[q[i] - 33] << std::endl;
    currP *= ( 1.0 - qualsMap[ q[i] - 33 ] );    
    encQ[i] = probToScoreChar(currP);
  }
  for (size_t i = k; i <= m-k; ++i) {
    currP /= (1.0 - qualsMap[q[i-k] - 33]);
    currP *= ( 1.0 - qualsMap[ q[i] - 33 ] );
    encQ[i] = probToScoreChar(currP);
    
  }
  for (size_t i = m - k + 1; i < m-1; ++i) {
    currP /= 1.0 - qualsMap[q[i-k] - 33];    
    encQ[i] = probToScoreChar(currP);
  }
  encQ[m-1] = q[m-1];
  return encQ;
  }*/

int main(int argc, char** argv) {

  std::cout << std::endl;
  
  // TODO: Add argument parsing code
  if (argc < 4) {
    std::cerr << "Not enough arguments" << std::endl;
    std::cerr << "  USAGE qComp [fastq_in] [fastq_out] [k]" << std::endl;
    exit(1);
  }

  init();
  
  std::string fastqFile(argv[1]);
  std::string outFilePath(argv[2]);
  size_t k = atoi(argv[3]);

  std::ifstream inFile(fastqFile.c_str());
  std::ofstream outFile(outFilePath.c_str());

  std::cout << "Input file               " << fastqFile << std::endl;
  std::cout << "Output file              " << outFilePath << std::endl;
  std::cout << "k                        " << k << std::endl;

  size_t readsCount = 0;
  Read r;
  while (!( r = nextRead(inFile) ).bases.empty()) {    
    r.qualities = reEncodeQuality(r.qualities, k);
    writeRead(r,outFile);
    readsCount++;
  }

  std::cout << "Reads count              " << readsCount << std::endl;

  inFile.close();
  outFile.close();

  std::cout << std::endl;

  // std::string QQ = "IIHHDI";
  // std::cout << QQ << std::endl;
  // std::cout << reEncodeQuality(QQ,1) << std::endl;


  
  return 0;
}
