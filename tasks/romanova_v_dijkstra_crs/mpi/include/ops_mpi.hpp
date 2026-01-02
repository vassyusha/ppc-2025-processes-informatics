#pragma once

#include <vector>

#include "romanova_v_dijkstra_crs/common/include/common.hpp"
#include "task/include/task.hpp"

namespace romanova_v_dijkstra_crs {

class RomanovaVDijkstraCrsMPI : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kMPI;
  }
  explicit RomanovaVDijkstraCrsMPI(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  void RemoveFromQueues(int vert);
  void CleanUpQueues();
  void UpdateQueues(int vert);

  Graph data_;
  std::vector<double> res_weights_;

  std::vector<double> min_in_;   // минимальный вес входящих ребер
  std::vector<double> min_out_;  // минимальный вес исходящих ребер

  int local_n_{};
  int st_vert_{};
  int en_vert_{};

  int delta_{};
  int extra_{};

  std::vector<double> local_d_;  // приоритеты для qin и qout будем вычислять динамически, используя min_in_, min_out_
  std::vector<bool> in_s_;
  std::vector<bool> visited_;

  MinHeap qd_;
  MinHeap qin_;
  MinHeap qout_;

  std::vector<bool> in_qd_;
  std::vector<bool> in_qin_;
  std::vector<bool> in_qout_;
};

}  // namespace romanova_v_dijkstra_crs
