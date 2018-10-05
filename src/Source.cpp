#include <stdlib.h>
#include <string>
#include <thread>
#include "MD.h"
#include "tinyxml2.h"

using namespace tinyxml2;

int main(int argc, char const* argv[]) {
  /* Path and name of the XML schema file */

  XMLDocument doc;
  doc.LoadFile("/home/gn/Code/MD-simulation/input_schema.xml");
}
