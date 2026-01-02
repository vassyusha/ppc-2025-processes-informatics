#pragma once

#include <cstddef>
#include <vector>

#include "romanova_v_dijkstra_crs/common/include/common.hpp"
#include "task/include/task.hpp"

namespace romanova_v_dijkstra_crs {

class RomanovaVDijkstraCrsSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit RomanovaVDijkstraCrsSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  Graph in_data_;
  std::vector<double> res_weights_;
};

}  // namespace romanova_v_dijkstra_crs
