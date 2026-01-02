#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <vector>

#include "romanova_v_dijkstra_crs/common/include/common.hpp"
#include "romanova_v_dijkstra_crs/mpi/include/ops_mpi.hpp"
#include "romanova_v_dijkstra_crs/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace romanova_v_dijkstra_crs {

class RomanovaVDijkstraCrsPerfTestProcesses : public ppc::util::BaseRunPerfTests<InType, OutType> {
  void SetUp() override {
    kVert_ = 1000;
    input_data_ = Graph();
    input_data_.vertices = kVert_;
    input_data_.offsets = std::vector<int>(input_data_.vertices + 1);
    input_data_.source = 0;

    for (int vert = 0; vert < kVert_; vert++) {
      for (int edg = 1; edg < 101; edg++) {
        int j = (vert + edg * 89) % kVert_;
        if (vert == j) {
          j = (j + 97) % kVert_;
        }
        input_data_.edges.push_back(j);
        auto weight = static_cast<double>(1 + ((vert + j + edg) % 37));
        input_data_.weights.push_back(weight);
        input_data_.offsets[vert + 1]++;
      }
    }

    for (int i = 0; i < kVert_; ++i) {
      input_data_.offsets[i + 1] += input_data_.offsets[i];
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != static_cast<size_t>(kVert_)) {
      return false;
    }

    if (abs(output_data[0]) > 1e-9) {
      return false;
    }

    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  int kVert_{};
};

TEST_P(RomanovaVDijkstraCrsPerfTestProcesses, Dijkstra) {
  ExecuteTest(GetParam());
}

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, RomanovaVDijkstraCrsMPI, RomanovaVDijkstraCrsSEQ>(
    PPC_SETTINGS_romanova_v_dijkstra_crs);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = RomanovaVDijkstraCrsPerfTestProcesses::CustomPerfTestName;
// NOLINTBEGIN(cppcoreguidelines-avoid-non-const-global-variables, modernize-type-traits, misc-use-anonymous-namespace)
INSTANTIATE_TEST_SUITE_P(RunModeTests, RomanovaVDijkstraCrsPerfTestProcesses, kGtestValues, kPerfTestName);
// NOLINTEND(cppcoreguidelines-avoid-non-const-global-variables, modernize-type-traits, misc-use-anonymous-namespace)
}  // namespace romanova_v_dijkstra_crs
