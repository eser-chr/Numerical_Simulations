#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "info.hpp"

/* All files : "info.cpp", "info.hpp" and "main.cpp" were laid out using
 * clang-format  */

void readFile(std::vector<Info> &vec);
void writeFile(std::vector<Info> &vec);
void simplify(std::vector<Info> &vec);

int main() {

  std::vector<Info> vec; // Initialising the vector who stores "Info" objects

  readFile(vec); // Reading the input files and storing information in "vec".

  

  simplify(vec); // Simplyfying the vector "vec".

  for (int i = 0; i < vec.size(); i++) {
    vec[i].print();
  }

  writeFile(vec); // Writing the output file from what is stored in "vec" after
                  // being simplified.

  return 0;
}

void readFile(std::vector<Info> &vec) {

  /*** Opening the files "transistors.dat and "specint.dat" ***/
  std::ifstream file1("transistors.dat");

  if (!file1) {
    std::cout << "Problem of opening file 1" << std::endl;
  }

  std::ifstream file2("specint.dat");

  if (!file2) {
    std::cout << "Problem of opening file 1" << std::endl;
  }

  /*** Initialisation of variables ***/
  double year_file1;
  double year_file2;

  double cpuCount;
  double perfNumber;

  std::string line1;
  std::string line2;

  /* Reading the first line from both files */
  std::getline(file1, line1);
  std::istringstream ligneStream(line1);
  ligneStream >> year_file1 >> cpuCount;

  std::getline(file2, line2);
  std::istringstream ligneStream2(line2);
  ligneStream2 >> year_file2 >> perfNumber;

  /* Reading the rest of the files */
  while (!(file1.eof() and file2.eof())) {

    /* To obtain a vector sorted in order of increasing date, the program must
     * compare the dates already read with each other and read from one or the
     * other file in order to maintain the correct order.*/

    if (year_file1 < year_file2) {
      // In this case -> reading from file 1

      if (line1[0] != '#') { // Ignoring line starting with a "#"
        Info tmp(year_file1, round(cpuCount * 1000), 0);        // Instianciating an objet with read informations
        vec.push_back(tmp); // Pushing this object in vec
      }

      /* Reading next line in file 1 */
      std::getline(file1, line1);
      std::istringstream ligneStream(line1);
      ligneStream >> year_file1 >> cpuCount;
    }

    else if (year_file1 > year_file2) {

      // In this case -> reading from file 2

      if (!(file2.eof())) { // Checking if file 2 is not over

        /* The following is the same as in the first case, but swapped with file
         * 2 */
        if (line2[0] != '#') {
          Info tmp(year_file2, 0, perfNumber / 1000);
          vec.push_back(tmp);
        }

        std::getline(file2, line2);
        std::istringstream ligneStream2(line2);
        ligneStream2 >> year_file2 >> perfNumber;
      }

      else {

        if (line1[0] != '#') {
          Info tmp(year_file1, round(cpuCount * 1000), 0);
          vec.push_back(tmp);
        }

        std::getline(file1, line1);
        std::istringstream ligneStream(line1);
        ligneStream >> year_file1 >> cpuCount;
      }
    }

    else if (year_file1 == year_file2) {

      /* In this case : reading in both files */

      if (line1[0] != '#' and line2[0] != '#') {
        Info tmp(year_file1, round(cpuCount * 1000), perfNumber / 1000);
        vec.push_back(tmp);
      }

      std::getline(file1, line1);
      std::istringstream ligneStream(line1);
      ligneStream >> year_file1 >> cpuCount;

      std::getline(file2, line2);
      std::istringstream ligneStream2(line2);
      ligneStream2 >> year_file2 >> perfNumber;
    }
  }
}

void writeFile(std::vector<Info> &vec) {

    /* Creating the ouput file as "output.txt" */
  std::ofstream output("output.txt");

    /* Checking if no problem happened */
  if (!(output.is_open())) {
    std::cout << "Problem during writing of the file" << std::endl;
  }

    /* For each object in "vec", writing its informations in the ouput file */
  size_t N = vec.size();

  for (int i = 0; i < N; i++) {
    output << vec[i].getYear() << ", " << vec[i].getCpuCount() << ", "
           << vec[i].getPerfNumber() << std::endl;
  }
  
}

void simplify(std::vector<Info> &vec) {

  size_t N = vec.size();

  for (int i = 1; i < N; i++) {

    int j = i;
    int count = 1;

    /* Counters for the number of informations counted in the average
     * calculation */
    int fac1 = 0;
    int fac2 = 0;

    /* Average values that will be computed before being stored */
    double avg_cpu_count = 0;
    double avg_perf_number = 0;

    while (vec[j - 1].getYear() == vec[j].getYear()) { // While the new year == last year --> count
      count++;
      j++;
    }

    // Computing the average values for the cpu count and the performance number
    for (int k = 0; k < count; k++) {

      if (vec[i + k - 1].getCpuCount() != 0) {
        avg_cpu_count += vec[i + k - 1].getCpuCount();
        fac1++;
      }
      if (vec[i + k - 1].getPerfNumber() != 0) {
        avg_perf_number += vec[i + k - 1].getPerfNumber();
        fac2++;
      }
    }

    // If there is more than one successive year
    if (count != 1) {

      /* Deletion of all elements of the vector which begins with the same year
        (except the first entry which will be modified) */
      vec.erase(vec.begin() + i, vec.begin() + i + count - 1);

      // Modification of the first element with the year
      if (fac1 != 0)
        vec[i - 1].setCpuCount(round(avg_cpu_count / fac1));
      if (fac2 != 0)
        vec[i - 1].setPerfNumber(avg_perf_number / fac2);

      // Modification of N and i so the loop does not get out of the size of the
      // vector
      N = N - count + 1;
      i = i - count;
    }
  }
}