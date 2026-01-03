#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <fstream>
#include <vector>

#include "romanova_v_dijkstra_crs/common/include/common.hpp"
#include "romanova_v_dijkstra_crs/mpi/include/ops_mpi.hpp"
#include "romanova_v_dijkstra_crs/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"
#include "util/include/util.hpp"

namespace romanova_v_dijkstra_crs {

class RomanovaVDijkstraCrsPerfTestProcesses : public ppc::util::BaseRunPerfTests<InType, OutType> {
  void SetUp() override {
    TestType params = "perf1000000";
    TestType abs_path = ppc::util::GetAbsoluteTaskPath(PPC_ID_romanova_v_dijkstra_crs, params);
    std::ifstream file(abs_path + ".txt");
    if (file.is_open()) {
      Graph gr;
      int edg = 0;
      file >> gr.vertices >> edg;
      gr.weights = std::vector<double>(edg);
      gr.edges = std::vector<int>(edg);
      gr.offsets = std::vector<int>(gr.vertices + 1);
      for (int i = 0; i < edg; i++) {
        file >> gr.weights[i];
      }
      for (int i = 0; i < edg; i++) {
        file >> gr.edges[i];
      }
      for (int i = 0; i < gr.vertices + 1; i++) {
        file >> gr.offsets[i];
      }
      file >> gr.source;
      file.close();
      input_data_ = gr;
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != static_cast<size_t>(1000000)) {
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
