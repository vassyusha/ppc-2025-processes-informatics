#pragma once

#include <mpi.h>

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

  void SetupMasterProcessData();
  void CalculateGlobalMinInOut();
  void BroadcastParameters();
  void SetupLocalVertexRange(int rank);
  void InitializeLocalArrays();
  void SetupMinInOutArrays(int rank, int n);
  std::vector<int> CalculateVertexSendCounts(int rank, int n);
  std::vector<int> CalculateVertexDisplacements(int rank, int n, const std::vector<int> &vert_sendcounts);
  void InitializeQueues();
  void BroadcastGraphData(int rank);

  void RecieveData(int &flag, MPI_Status &status);
  void MakeLocalR(std::vector<int> &local_r, double global_l, double global_m);
  void ProcessLocalR(std::vector<int> &local_r, std::vector<MPI_Request> &dist_requests,
                     std::vector<MPI_Request> &vertex_requests, int &flag, MPI_Status &status);
  bool IsGlobalStop();
  double GetGlobalMin(MinHeap &q);

  Graph data_;
  std::vector<double> res_weights_;

  std::vector<double> glob_min_in_;
  std::vector<double> glob_min_out_;

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
