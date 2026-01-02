#include "romanova_v_dijkstra_crs/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <queue>
#include <tuple>
#include <vector>

#include "romanova_v_dijkstra_crs/common/include/common.hpp"

namespace romanova_v_dijkstra_crs {

RomanovaVDijkstraCrsSEQ::RomanovaVDijkstraCrsSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput();
}

bool RomanovaVDijkstraCrsSEQ::ValidationImpl() {
  if (GetInput().vertices <= 0) {
    return false;
  }
  if (GetInput().offsets.size() - 1 != GetInput().vertices) {
    return false;
  }
  if (GetInput().source < 0 || GetInput().source >= GetInput().vertices) {
    return false;
  }
  for (size_t i = 0; i < GetInput().weights.size(); i++) {
    if (GetInput().weights[i] < 0) {
      return false;
    }
  }
  return true;
}

bool RomanovaVDijkstraCrsSEQ::PreProcessingImpl() {
  in_data_ = GetInput();
  return true;
}

bool RomanovaVDijkstraCrsSEQ::RunImpl() {
  int n = in_data_.vertices;
  int source = in_data_.source;
  res_weights_ = std::vector<double>(n, std::numeric_limits<double>::infinity());
  res_weights_[source] = 0.0;

  std::vector<bool> visited(n, false);

  std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<>> pq;
  pq.emplace(std::make_pair(0.0, source));

  while (!pq.empty()) {
    auto [curr_dist, vert] = pq.top();
    pq.pop();

    if (visited[vert]) {
      continue;
    }
    visited[vert] = true;

    int st = in_data_.offsets[vert];
    int end = in_data_.offsets[vert + 1];
    for (int i = st; i < end; i++) {
      if (visited[in_data_.edges[i]]) {
        continue;
      }
      if (curr_dist + in_data_.weights[i] < res_weights_[in_data_.edges[i]]) {
        res_weights_[in_data_.edges[i]] = curr_dist + in_data_.weights[i];
        pq.emplace(std::make_pair(res_weights_[in_data_.edges[i]], in_data_.edges[i]));
      }
    }
  }

  for (int i = 0; i < in_data_.vertices; i++) {
    if (res_weights_[i] == std::numeric_limits<double>::infinity()) {
      res_weights_[i] = -1;
    }
  }

  return true;
}

bool RomanovaVDijkstraCrsSEQ::PostProcessingImpl() {
  GetOutput() = res_weights_;
  return true;
}

}  // namespace romanova_v_dijkstra_crs
