class Info {
public:
  // Constructor
  Info();
  Info(long year, double cpuCount, double perfNumber);

  // Getter fonctions
  long getYear() const;
  double getCpuCount() const;
  double getPerfNumber() const;

  // Setter fonctions 
  void setYear(long newYear);
  void setCpuCount(double newCpuCount);
  void setPerfNumber(double newPerfNumber);

  // Others fonctions
  void modify(long year, double cpuCount, double perfNumber);
  void print();

private:
  long _year;
  double _cpuCount;
  double _perfNumber;
};