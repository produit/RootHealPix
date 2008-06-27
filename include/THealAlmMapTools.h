#ifndef T_HEAL_AML_MAP_TOOLS_H
#define T_HEAL_AML_MAP_TOOLS_H

#include <vector>

#include "THealAlm.h"

class THealPix;

namespace THealAlmMapTools {
  template<typename T> void Map2Alm (const THealPix &map, THealAlm<T> &alm, const std::vector<double> &weight, bool add_alm);
  template<typename T> void Alm2Map (const THealAlm<T> &alm, THealPix &map);
}

#endif // T_HEAL_AML_MAP_TOOLS_H
