#include "romanova_v_dijkstra_crs/mpi/include/ops_mpi.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

#include "romanova_v_dijkstra_crs/common/include/common.hpp"

namespace romanova_v_dijkstra_crs {

void RomanovaVDijkstraCrsMPI::RemoveFromQueues(int vert) {
  in_qd_[vert] = false;
  in_qin_[vert] = false;
  in_qout_[vert] = false;
}

void RomanovaVDijkstraCrsMPI::CleanUpQueues() {
  while (!qd_.empty()) {
    auto item = qd_.top();
    if (!in_s_[item.second]) {
      break;
    }
    qd_.pop();
    in_qd_[item.second] = false;
  }

  while (!qin_.empty()) {
    auto item = qin_.top();
    if (!in_s_[item.second]) {
      break;
    }
    qin_.pop();
    in_qin_[item.second] = false;
  }

  while (!qout_.empty()) {
    auto item = qout_.top();
    if (!in_s_[item.second]) {
      break;
    }
    qout_.pop();
    in_qout_[item.second] = false;
  }
}

void RomanovaVDijkstraCrsMPI::UpdateQueues(int vert) {
  qd_.emplace(local_d_[vert], vert);
  in_qd_[vert] = true;

  qin_.emplace(local_d_[vert] - min_in_[vert], vert);
  in_qin_[vert] = true;

  qout_.emplace(local_d_[vert] + min_out_[vert], vert);
  in_qout_[vert] = true;
}

void RomanovaVDijkstraCrsMPI::RecieveData(int &flag, MPI_Status &status) {
  MPI_Iprobe(MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &flag, &status);

  while (flag != 0) {
    struct SendData {
      double distance;
      int vertex;
    } recieved_data{};

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Recv(&recieved_data, sizeof(SendData), MPI_BYTE, status.MPI_SOURCE, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    double new_dist = recieved_data.distance;
    int glob_v = recieved_data.vertex;
    if (st_vert_ <= glob_v && glob_v < en_vert_) {
      if (new_dist < local_d_[glob_v - st_vert_]) {
        local_d_[glob_v - st_vert_] = new_dist;
        UpdateQueues(glob_v - st_vert_);
      }
    }

    MPI_Iprobe(MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &flag, &status);
  }
}

void RomanovaVDijkstraCrsMPI::WaitRequests(std::vector<MPI_Request> &send_requests) {
  if (!send_requests.empty()) {
    MPI_Waitall(static_cast<int>(send_requests.size()), send_requests.data(), MPI_STATUSES_IGNORE);
  }
}

void RomanovaVDijkstraCrsMPI::MakeLocalR(std::vector<int> &local_r, double global_l, double global_m) {
  for (int i = 0; i < local_n_; i++) {
    if (!in_s_[i]) {
      bool cond1 = (local_d_[i] <= global_l + 1e-9);
      bool cond2 = (local_d_[i] - min_in_[i] <= global_m + 1e-9);

      if (cond1 || cond2) {
        local_r.push_back(i);
        in_s_[i] = true;
      }
    }
  }
}

void RomanovaVDijkstraCrsMPI::ProcessLocalR(std::vector<int> &local_r, int &flag, MPI_Status &status) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::vector<MPI_Request> send_requests;
  for (int u : local_r) {
    int start = data_.offsets[u + st_vert_];
    int end = data_.offsets[u + st_vert_ + 1];

    // NOLINTBEGIN(clang-analyzer-optin.mpi.MPI-Checker)
    for (int j = start; j < end; j++) {
      int glob_v = data_.edges[j];
      double weight = data_.weights[j];

      double new_dist = local_d_[u] + weight;

      int owner = (glob_v < extra_ * (delta_ + 1) ? glob_v / (delta_ + 1)
                                                  : extra_ + ((glob_v - (delta_ + 1) * extra_) / delta_));

      if (owner == rank) {
        if (new_dist < local_d_[glob_v - st_vert_]) {
          local_d_[glob_v - st_vert_] = new_dist;
          UpdateQueues(glob_v - st_vert_);
        }
      } else {
        struct SendData {
          double distance;
          int vertex;
        } send_data{.distance = new_dist, .vertex = glob_v};

        MPI_Request req = MPI_REQUEST_NULL;
        MPI_Isend(&send_data, sizeof(SendData), MPI_BYTE, owner, 2, MPI_COMM_WORLD, &req);
        send_requests.push_back(req);
      }
    }
    // NOLINTEND(clang-analyzer-optin.mpi.MPI-Checker)
    RecieveData(flag, status);
  }
  WaitRequests(send_requests);
}

std::vector<int> RomanovaVDijkstraCrsMPI::IsGlobalStop() {
  int local_has_work = (!qd_.empty() || !qin_.empty() || !qout_.empty()) ? 1 : 0;

  int has_pending_msgs = 0;
  MPI_Status temp_status;
  MPI_Iprobe(MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &has_pending_msgs, &temp_status);

  while (has_pending_msgs != 0) {
    RecieveData(has_pending_msgs, temp_status);
    local_has_work = (!qd_.empty() || !qin_.empty() || !qout_.empty()) ? 1 : 0;
    MPI_Iprobe(MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &has_pending_msgs, &temp_status);
  }

  std::vector<int> status = {local_has_work, has_pending_msgs};
  return status;
}

double RomanovaVDijkstraCrsMPI::GetGlobalMin(MinHeap &q) {
  double local = std::numeric_limits<double>::infinity();
  if (!q.empty()) {
    local = q.top().first;
  }

  double global = NAN;
  MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  return global;
}

void RomanovaVDijkstraCrsMPI::SetupMasterProcessData() {
  data_ = GetInput();

  int n = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &n);

  delta_ = data_.vertices / n;
  extra_ = data_.vertices % n;

  CalculateGlobalMinInOut();
}

void RomanovaVDijkstraCrsMPI::CalculateGlobalMinInOut() {
  int vertices = data_.vertices;
  std::vector<double> glob_min_in(vertices, std::numeric_limits<double>::infinity());
  std::vector<double> glob_min_out(vertices, std::numeric_limits<double>::infinity());

  for (size_t i = 0; i < data_.offsets.size() - 1; i++) {
    for (int j = data_.offsets[i]; j < data_.offsets[i + 1]; j++) {
      glob_min_out[i] = std::min(glob_min_out[i], data_.weights[j]);

      int target_vertex = data_.edges[j];
      glob_min_in[target_vertex] = std::min(glob_min_in[target_vertex], data_.weights[j]);
    }
  }

  glob_min_in_ = std::move(glob_min_in);
  glob_min_out_ = std::move(glob_min_out);
}

void RomanovaVDijkstraCrsMPI::BroadcastParameters() {
  MPI_Bcast(&delta_, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&extra_, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

void RomanovaVDijkstraCrsMPI::SetupLocalVertexRange(int rank) {
  local_n_ = delta_ + (rank < extra_ ? 1 : 0);
  st_vert_ = (rank * delta_) + (rank < extra_ ? rank : extra_);
  en_vert_ = st_vert_ + local_n_;
}

void RomanovaVDijkstraCrsMPI::InitializeLocalArrays() {
  local_d_.assign(local_n_, std::numeric_limits<double>::infinity());
  in_s_.assign(local_n_, false);
  visited_.assign(local_n_, false);

  min_in_.assign(local_n_, 0.0);
  min_out_.assign(local_n_, 0.0);
}

void RomanovaVDijkstraCrsMPI::SetupMinInOutArrays(int rank, int n) {
  std::vector<int> vert_sendcounts = CalculateVertexSendCounts(rank, n);
  std::vector<int> vert_displs = CalculateVertexDisplacements(rank, n, vert_sendcounts);

  MPI_Scatterv(glob_min_in_.data(), vert_sendcounts.data(), vert_displs.data(), MPI_DOUBLE, min_in_.data(), local_n_,
               MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(glob_min_out_.data(), vert_sendcounts.data(), vert_displs.data(), MPI_DOUBLE, min_out_.data(), local_n_,
               MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

std::vector<int> RomanovaVDijkstraCrsMPI::CalculateVertexSendCounts(int rank, int n) const {
  std::vector<int> vert_sendcounts(n, delta_);

  if (rank == 0) {
    vert_sendcounts[0] += (extra_ > 0 ? 1 : 0);
    for (int i = 1; i < n; i++) {
      vert_sendcounts[i] += (i < extra_ ? 1 : 0);
    }
  }

  return vert_sendcounts;
}

std::vector<int> RomanovaVDijkstraCrsMPI::CalculateVertexDisplacements(int rank, int n,
                                                                       const std::vector<int> &vert_sendcounts) {
  std::vector<int> vert_displs(n, 0);

  if (rank == 0) {
    for (int i = 1; i < n; i++) {
      vert_displs[i] = vert_displs[i - 1] + vert_sendcounts[i - 1];
    }
  }

  return vert_displs;
}

void RomanovaVDijkstraCrsMPI::InitializeQueues() {
  in_qd_.assign(local_n_, false);
  in_qin_.assign(local_n_, false);
  in_qout_.assign(local_n_, false);

  qd_ = MinHeap();
  qin_ = MinHeap();
  qout_ = MinHeap();
}

void RomanovaVDijkstraCrsMPI::BroadcastGraphData(int rank) {
  int edg_sz = static_cast<int>(data_.edges.size());

  MPI_Bcast(&data_.vertices, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&data_.source, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&edg_sz, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank != 0) {
    data_.edges.assign(edg_sz, 0);
    data_.weights.assign(edg_sz, 0.0);
    data_.offsets.assign(data_.vertices + 1, 0);
  }

  MPI_Bcast(data_.edges.data(), static_cast<int>(data_.edges.size()), MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(data_.weights.data(), static_cast<int>(data_.weights.size()), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(data_.offsets.data(), static_cast<int>(data_.offsets.size()), MPI_INT, 0, MPI_COMM_WORLD);
}

RomanovaVDijkstraCrsMPI::RomanovaVDijkstraCrsMPI(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput();
}

bool RomanovaVDijkstraCrsMPI::ValidationImpl() {
  bool status = true;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    status = status && (GetInput().vertices > 0);
    status = status && (GetInput().offsets.size() - 1 == static_cast<size_t>(GetInput().vertices));
    status = status && (GetInput().source >= 0 && GetInput().source < GetInput().vertices);
    for (double weight : GetInput().weights) {
      status = status && (weight >= 1e-9);
    }
  }

  MPI_Bcast(&status, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
  return status;
}

bool RomanovaVDijkstraCrsMPI::PreProcessingImpl() {
  delta_ = 0;
  extra_ = 0;

  int rank = 0;
  int n = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n);

  if (rank == 0) {
    SetupMasterProcessData();
  }

  BroadcastParameters();
  SetupLocalVertexRange(rank);
  InitializeLocalArrays();
  SetupMinInOutArrays(rank, n);
  InitializeQueues();
  BroadcastGraphData(rank);

  return true;
}

void RomanovaVDijkstraCrsMPI::InitializeSource() {
  if (st_vert_ <= data_.source && data_.source < en_vert_) {
    local_d_[data_.source - st_vert_] = 0.0;
    UpdateQueues(data_.source - st_vert_);
  }
}

bool RomanovaVDijkstraCrsMPI::RunImpl() {
  int rank = 0;
  int n = 0;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n);

  InitializeSource();

  bool global_stop = false;

  while (!global_stop) {
    MPI_Status status;
    int flag = 0;

    RecieveData(flag, status);

    double global_l = GetGlobalMin(qout_);
    double global_m = GetGlobalMin(qd_);

    std::vector<int> local_r;
    MakeLocalR(local_r, global_l, global_m);

    for (int v : local_r) {
      RemoveFromQueues(v);
    }

    ProcessLocalR(local_r, flag, status);

    CleanUpQueues();
    MPI_Barrier(MPI_COMM_WORLD);

    std::vector<int> stop = IsGlobalStop();
    int global_has_work = 0;
    MPI_Allreduce(&stop.data(), &global_has_work, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);

    int global_has_pending = 0;
    MPI_Allreduce(&stop[1], &global_has_pending, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);

    MPI_Status temp_status;
    if (global_has_work == 0 && global_has_pending == 0) {
      MPI_Iprobe(MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &stop[1], &temp_status);
      global_stop = (stop[1] == 0);
    }

    // global_stop = IsGlobalStop();
  }

  MPI_Status final_status;
  int final_flag = 0;
  MPI_Iprobe(MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &final_flag, MPI_STATUS_IGNORE);
  if (final_flag != 0) {
    RecieveData(final_flag, final_status);
  }

  std::vector<int> vert_sendcounts(n, delta_);
  std::vector<int> vert_displs(n, 0);

  if (rank == 0) {
    vert_sendcounts = CalculateVertexSendCounts(rank, n);
    vert_displs = CalculateVertexDisplacements(rank, n, vert_sendcounts);
  }

  res_weights_.assign(data_.vertices, 0.0);

  MPI_Gatherv(local_d_.data(), local_n_, MPI_DOUBLE, res_weights_.data(), vert_sendcounts.data(), vert_displs.data(),
              MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (rank == 0) {
    for (int i = 0; i < data_.vertices; i++) {
      if (res_weights_[i] == std::numeric_limits<double>::infinity()) {
        res_weights_[i] = -1;
      }
    }
  }
  MPI_Bcast(res_weights_.data(), data_.vertices, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  return true;
}

bool RomanovaVDijkstraCrsMPI::PostProcessingImpl() {
  GetOutput() = res_weights_;
  return true;
}

}  // namespace romanova_v_dijkstra_crs
