#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#include "romanova_v_dijkstra_crs/common/include/common.hpp"
#include "romanova_v_dijkstra_crs/mpi/include/ops_mpi.hpp"
#include "romanova_v_dijkstra_crs/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace romanova_v_dijkstra_crs {

class RomanovaVDijkstraCrsFuncTestsProcesses : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return test_param;
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    std::string abs_path = ppc::util::GetAbsoluteTaskPath(PPC_ID_romanova_v_dijkstra_crs, params);
    std::ifstream file(abs_path + ".txt");
    if (file.is_open()) {
      input_data_ = Graph();
      int edg = 0;
      file >> input_data_.vertices >> edg;
      input_data_.weights = std::vector<double>(edg);
      input_data_.edges = std::vector<int>(edg);
      input_data_.offsets = std::vector<int>(input_data_.vertices + 1);
      for (int i = 0; i < edg; i++) {
        file >> input_data_.weights[i];
      }
      for (int i = 0; i < edg; i++) {
        file >> input_data_.edges[i];
      }
      for (int i = 0; i < input_data_.vertices + 1; i++) {
        file >> input_data_.offsets[i];
      }
      file >> input_data_.source;

      exp_answer_ = OutType(input_data_.vertices);
      for (int i = 0; i < input_data_.vertices; i++) {
        file >> exp_answer_[i];
      }
      file.close();
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    for (size_t i = 0; i < exp_answer_.size(); i++) {
      if (abs(output_data[i] - exp_answer_[i]) > 1e-9) {
        return false;
      }
    }

    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType exp_answer_;
};

namespace {

TEST_P(RomanovaVDijkstraCrsFuncTestsProcesses, Dijkstra) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 9> kTestParam = {"trivialTest",     "smallTest",       "simpleTest",
                                            "disconGraphTest", "linGraphTest",    "cycleGraphTest",
                                            "complGraphTest",  "severalWaysTest", "longerWayWithLessCostTest"};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<RomanovaVDijkstraCrsMPI, InType>(kTestParam, PPC_SETTINGS_romanova_v_dijkstra_crs),
    ppc::util::AddFuncTask<RomanovaVDijkstraCrsSEQ, InType>(kTestParam, PPC_SETTINGS_romanova_v_dijkstra_crs));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName =
    RomanovaVDijkstraCrsFuncTestsProcesses::PrintFuncTestName<RomanovaVDijkstraCrsFuncTestsProcesses>;
// NOLINTBEGIN(cppcoreguidelines-avoid-non-const-global-variables, modernize-type-traits)
INSTANTIATE_TEST_SUITE_P(Tests, RomanovaVDijkstraCrsFuncTestsProcesses, kGtestValues, kPerfTestName);
// NOLINTEND(cppcoreguidelines-avoid-non-const-global-variables, modernize-type-traits)
}  // namespace

}  // namespace romanova_v_dijkstra_crs
