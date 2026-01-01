#pragma once

#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace romanova_v_dijkstra_crs {

struct Graph{
    std::vector<double> weights;
    std::vector<int> edges;
    std::vector<int> offsets;

    int vertices;
    int source;
};

using InType = Graph;
using OutType = std::vector<double>;
using TestType = std::string;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace romanova_v_dijkstra_crs
