#include "romanova_v_dijkstra_crs/mpi/include/ops_mpi.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <tuple>
#include <utility>
#include <vector>

#include "romanova_v_dijkstra_crs/common/include/common.hpp"

namespace romanova_v_dijkstra_crs {

RomanovaVDijkstraCrsMPI::RomanovaVDijkstraCrsMPI(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput();
}

bool RomanovaVDijkstraCrsMPI::ValidationImpl() {
  return true;
}

bool RomanovaVDijkstraCrsMPI::PreProcessingImpl() {
  return true;
}

bool RomanovaVDijkstraCrsMPI::RunImpl() {
  return true;
}

bool RomanovaVDijkstraCrsMPI::PostProcessingImpl() {
  return true;
}

}  // namespace romanova_v_dijkstra_crs
