#include "info.hpp"
#include <iostream>

/********* Constructors *********/

Info::Info() {
  _year = 0;
  _cpuCount = 0;
  _perfNumber = 0;
}

Info::Info(long year, double cpuCount, double perfNumber) {
  _year = year;
  _cpuCount = cpuCount;
  _perfNumber = perfNumber;
}

/********* Getter fonctions *********/

long Info::getYear() const { return _year; }

double Info::getCpuCount() const { return _cpuCount; }

double Info::getPerfNumber() const { return _perfNumber; }

/********* Setter fonctions *********/

void Info::setYear(long newYear) { _year = newYear; }

void Info::setCpuCount(double newCpuCount) { _cpuCount = newCpuCount; }

void Info::setPerfNumber(double newPerfNumber) { _perfNumber = newPerfNumber; }

void Info::modify(long year, double cpuCount, double perfNumber) {
  _year = year;
  _cpuCount = cpuCount;
  _perfNumber = perfNumber;
}

/********* Other fonctions *********/

void Info::print() {
  std::cout << _year << " / " << _cpuCount << " / " << _perfNumber << std::endl;
}