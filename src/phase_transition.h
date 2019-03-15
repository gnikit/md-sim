#include "MD.h"

class phase_transition : public MD {
 private:
  /* data */
 public:
    using MD::MD;

  void crystallisation(double DENSITY, double FINAL_DENSITY, double DENSITY_INC,
                       double TEMPERATURE, double POWER, double A_CST, std::string pp_type);
};