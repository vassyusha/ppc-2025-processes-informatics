#include "romanova_v_dijkstra_crs/mpi/include/ops_mpi.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <queue>
#include <tuple>
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
    if (in_qd_[item.second]) {
      break;
    }
    qd_.pop();
  }

  while (!qin_.empty()) {
    auto item = qin_.top();
    if (in_qin_[item.second]) {
      break;
    }
    qin_.pop();
  }

  while (!qout_.empty()) {
    auto item = qout_.top();
    if (in_qout_[item.second]) {
      break;
    }
    qout_.pop();
  }
}

void RomanovaVDijkstraCrsMPI::UpdateQueues(int vert) {
  if (!in_qd_[vert]) {
    qd_.emplace(local_d_[vert], vert);
    in_qd_[vert] = true;
  }

  if (!in_qin_[vert]) {
    qin_.emplace(local_d_[vert] - min_in_[vert], vert);
    in_qin_[vert] = true;
  }

  if (!in_qout_[vert]) {
    qout_.emplace(local_d_[vert] + min_out_[vert], vert);
    in_qout_[vert] = true;
  }
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
    status = status && (GetInput().offsets.size() - 1 == GetInput().vertices);
    status = status && (GetInput().source >= 0 && GetInput().source < GetInput().vertices);
    for (size_t i = 0; i < GetInput().weights.size(); i++) {
      status = status && (GetInput().weights[i] >= 1e-9);
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

  std::vector<double> glob_min_in;
  std::vector<double> glob_min_out;

  data_ = Graph();

  if (rank == 0) {
    data_ = GetInput();

    delta_ = data_.vertices / n;
    extra_ = data_.vertices % n;

    glob_min_in = std::vector<double>(data_.vertices, std::numeric_limits<double>::infinity());
    glob_min_out = std::vector<double>(data_.vertices, std::numeric_limits<double>::infinity());

    for (size_t i = 0; i < data_.offsets.size() - 1; i++) {
      for (size_t j = data_.offsets[i]; j < data_.offsets[i + 1]; j++) {
        if (data_.weights[j] < glob_min_out[i]) {
          glob_min_out[i] = data_.weights[j];
        }
        if (data_.weights[j] < glob_min_in[data_.edges[j]]) {
          glob_min_in[data_.edges[j]] = data_.weights[j];
        }
      }
    }
  }

  MPI_Bcast(&delta_, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&extra_, 1, MPI_INT, 0, MPI_COMM_WORLD);

  local_n_ = delta_ + (rank < extra_ ? 1 : 0);
  st_vert_ = rank * delta_ + (rank < extra_ ? rank : extra_);
  en_vert_ = st_vert_ + local_n_;

  local_d_.assign(local_n_, std::numeric_limits<double>::infinity());
  in_s_.assign(local_n_, false);
  visited_.assign(local_n_, false);

  min_in_.assign(local_n_, 0.0);
  min_out_.assign(local_n_, 0.0);

  std::vector<int> vert_sendcounts;
  std::vector<int> vert_displs;

  if (rank == 0) {
    vert_sendcounts.assign(n, delta_);
    vert_sendcounts[0] += (extra_ > 0 ? 1 : 0);
    vert_displs.assign(n, 0);

    for (int i = 1; i < n; i++) {
      vert_sendcounts[i] += (i < extra_ ? 1 : 0);
      vert_displs[i] = vert_displs[i - 1] + vert_sendcounts[i - 1];
    }
  }

  MPI_Scatterv(glob_min_in.data(), vert_sendcounts.data(), vert_displs.data(), MPI_DOUBLE, min_in_.data(), local_n_,
               MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(glob_min_out.data(), vert_sendcounts.data(), vert_displs.data(), MPI_DOUBLE, min_out_.data(), local_n_,
               MPI_DOUBLE, 0, MPI_COMM_WORLD);

  in_qd_.assign(local_n_, false);
  in_qin_.assign(local_n_, false);
  in_qout_.assign(local_n_, false);

  qd_ = MinHeap();
  qin_ = MinHeap();
  qout_ = MinHeap();

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

  return true;
}

bool RomanovaVDijkstraCrsMPI::RunImpl() {
  int rank = 0;
  int n = 0;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n);

  if (st_vert_ <= data_.source && data_.source < en_vert_) {
    local_d_[data_.source - st_vert_] = 0.0;
    UpdateQueues(data_.source - st_vert_);
  }

  bool global_stop = false;

  while (!global_stop) {
    MPI_Status status;
    int flag;
    MPI_Iprobe(MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &flag, &status);

    while (flag) {
      double new_dist = 0.0;
      int glob_v = 0;
      MPI_Recv(&new_dist, 1, MPI_DOUBLE, status.MPI_SOURCE, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&glob_v, 1, MPI_INT, status.MPI_SOURCE, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      if (st_vert_ <= glob_v && glob_v < en_vert_) {
        if (new_dist < local_d_[glob_v - st_vert_]) {
          local_d_[glob_v - st_vert_] = new_dist;
          UpdateQueues(glob_v - st_vert_);
        }
      }

      MPI_Iprobe(MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &flag, &status);
    }

    double local_l = std::numeric_limits<double>::infinity();
    if (!qout_.empty()) {
      local_l = qout_.top().first;
    }

    double local_m = std::numeric_limits<double>::infinity();
    if (!qd_.empty()) {
      local_m = qd_.top().first;
    }

    double global_l;
    double global_m;
    MPI_Allreduce(&local_l, &global_l, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&local_m, &global_m, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    std::vector<int> local_r;
    for (int i = 0; i < local_n_; i++) {
      if (!visited_[i] && !in_s_[i]) {
        bool cond1 = (local_d_[i] <= global_l);
        bool cond2 = (local_d_[i] - min_in_[i] <= global_m);

        if (cond1 || cond2) {
          local_r.push_back(i);
          in_s_[i] = true;
          visited_[i] = true;
        }
      }
    }

    for (int v : local_r) {
      RemoveFromQueues(v);
    }

    std::vector<MPI_Request> send_requests;

    for (int u : local_r) {
      size_t start = static_cast<size_t>(data_.offsets[u + st_vert_]);
      size_t end = static_cast<size_t>(data_.offsets[u + st_vert_ + 1]);
      for (size_t j = start; j < end; j++) {
        int glob_v = data_.edges[j];
        double weight = data_.weights[j];

        double new_dist = local_d_[u] + weight;

        int owner = (glob_v < extra_ * (delta_ + 1) ? glob_v / (delta_ + 1)
                                                    : extra_ + (glob_v - (delta_ + 1) * extra_) / delta_);
        if (owner == rank) {
          if (new_dist < local_d_[glob_v - st_vert_]) {
            local_d_[glob_v - st_vert_] = new_dist;
            UpdateQueues(glob_v - st_vert_);
          }
        } else {
          MPI_Request req1, req2;
          MPI_Isend(&new_dist, 1, MPI_DOUBLE, owner, 2, MPI_COMM_WORLD, &req1);
          MPI_Isend(&glob_v, 1, MPI_INT, owner, 3, MPI_COMM_WORLD, &req2);
          send_requests.push_back(req1);
          send_requests.push_back(req2);
        }
      }

      MPI_Iprobe(MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &flag, &status);

      while (flag) {
        double new_dist = 0.0;
        int glob_v = 0;
        MPI_Recv(&new_dist, 1, MPI_DOUBLE, status.MPI_SOURCE, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&glob_v, 1, MPI_INT, status.MPI_SOURCE, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        if (new_dist < local_d_[glob_v - st_vert_]) {
          local_d_[glob_v - st_vert_] = new_dist;
          UpdateQueues(glob_v - st_vert_);
        }

        MPI_Iprobe(MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &flag, &status);
      }
    }

    if (!send_requests.empty()) {
      MPI_Waitall(static_cast<int>(send_requests.size()), send_requests.data(), MPI_STATUSES_IGNORE);
    }
    CleanUpQueues();
    MPI_Barrier(MPI_COMM_WORLD);

    int local_has_work = (!qd_.empty()) ? 1 : 0;
    int global_has_work;
    MPI_Allreduce(&local_has_work, &global_has_work, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);

    global_stop = (global_has_work == 0);
  }

  MPI_Status final_status;
  int final_flag;

  MPI_Iprobe(MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &final_flag, &final_status);
  while (final_flag) {
    double new_dist = 0.0;
    int glob_v = 0;
    MPI_Recv(&new_dist, 1, MPI_DOUBLE, final_status.MPI_SOURCE, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&glob_v, 1, MPI_INT, final_status.MPI_SOURCE, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if (st_vert_ <= glob_v && glob_v < en_vert_) {
      if (new_dist < local_d_[glob_v - st_vert_]) {
        local_d_[glob_v - st_vert_] = new_dist;
      }
    }

    MPI_Iprobe(MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &final_flag, &final_status);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  std::vector<int> vert_sendcounts;
  std::vector<int> vert_displs;

  if (rank == 0) {
    vert_sendcounts.assign(n, delta_);
    vert_sendcounts[0] += (extra_ > 0 ? 1 : 0);
    vert_displs.assign(n, 0);

    for (int i = 1; i < n; i++) {
      vert_sendcounts[i] += (i < extra_ ? 1 : 0);
      vert_displs[i] = vert_displs[i - 1] + vert_sendcounts[i - 1];
    }
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
