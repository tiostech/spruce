/*-------------------------------------------------------------------------------
 This file is part of spruce.

 Copyright (c) [2014-2018] [Marvin N. Wright]

 This software may be modified and distributed under the terms of the MIT license.

 Please note that the C++ core of spruce is distributed under MIT license and the
 R package "spruce" under GPL3 license.
 #-------------------------------------------------------------------------------*/

#ifndef TREEREGRESSION_H_
#define TREEREGRESSION_H_

#include <vector>

#include "globals.h"
#include "Tree.h"

namespace spruce
{

  class TreeRegression : public Tree
  {
  public:
    TreeRegression() = default;

    // Create from loaded forest
    TreeRegression(std::vector<std::vector<size_t>> &child_nodeIDs, std::vector<size_t> &split_varIDs, std::vector<double> &split_values);

    TreeRegression(const TreeRegression &) = delete;
    TreeRegression &operator=(const TreeRegression &) = delete;

    virtual ~TreeRegression() override = default;

    void allocateMemory() override;

    double estimate(size_t nodeID);
    void computePermutationImportanceInternal(std::vector<std::vector<size_t>> *permutations);
    void appendToFileInternal(std::ofstream &file) override;

    double getPrediction(size_t sampleID) const
    {
      size_t terminal_nodeID = prediction_terminal_nodeIDs[sampleID];
      return (split_values[terminal_nodeID]); // git fitted value 
    }

    size_t getPredictionTerminalNodeID(size_t sampleID) const
    {
      return prediction_terminal_nodeIDs[sampleID];
    }
    
    // .......
    double crystal_cost_core(double z, double y, double estimate);
    double crystal_cost(double z, double y, double estimate);
    double crystal_cost(std::vector<double> z_values, std::vector<double> y_values, double estimate);
    double crystal_fit(size_t nodeID);
    double crystal_fit(std::vector<double> z_values, std::vector<double> y_values);
    double crystal_fast_fit(std::vector<double> z_values, std::vector<double> y_values);
    double crystal_fast_fit(size_t nodeID); 
    
    bool findBestSplitCrystal(size_t nodeID, std::vector<size_t> &possible_split_varIDs);
    void findBestSplitValueCrystal(size_t nodeID, size_t varID, size_t num_samples_node, double &best_value, size_t &best_varID, double &best_decrease);
    void findBestSplitValueCrystal(size_t nodeID, size_t varID, size_t num_samples_node, double &best_value, size_t &best_varID, double &best_decrease, std::vector<double> possible_split_values, std::vector<double> &sums_right, std::vector<size_t> &n_right);
    bool splitNodeInternalCrystal(size_t nodeID, std::vector<size_t> &possible_split_varIDs) override;
    //........
    
    double absolute_cost(std::vector<double> z_values, std::vector<double> y_values, double estimate);
    double absolute_cost_fit(std::vector<double> z_values, std::vector<double> y_values);
    double absolute_cost_fit(size_t nodeID);

    bool findBestSplitAbsoluteCost(size_t nodeID, std::vector<size_t> &possible_split_varIDs);
    void findBestSplitValueAbsoluteCost(size_t nodeID, size_t varID, size_t num_samples_node, double &best_value, size_t &best_varID, double &best_decrease);
    void findBestSplitValueAbsoluteCost(size_t nodeID, size_t varID, size_t num_samples_node, double &best_value, size_t &best_varID, double &best_decrease, std::vector<double> possible_split_values, std::vector<double> &sums_right, std::vector<size_t> &n_right);
    bool splitNodeInternalAbsoluteCost(size_t nodeID, std::vector<size_t> &possible_split_varIDs) override;
    
    bool splitNodeInternalV2(size_t nodeID, std::vector<size_t> &possible_split_varIDs);
    
    
  // private:
    bool splitNodeInternal(size_t nodeID, std::vector<size_t> &possible_split_varIDs) override;
    void createEmptyNodeInternal() override;

    double computePredictionAccuracyInternal(std::vector<double> *prediction_error_casewise) override;

    // Called by splitNodeInternal(). Sets split_varIDs and split_values.



    bool findBestSplit(size_t nodeID, std::vector<size_t> &possible_split_varIDs);
    void findBestSplitValueSmallQ(size_t nodeID, size_t varID, double sum_node, size_t num_samples_node, double &best_value, size_t &best_varID, double &best_decrease);
    void findBestSplitValueSmallQ(size_t nodeID, size_t varID, double sum_node, size_t num_samples_node, double &best_value, size_t &best_varID, double &best_decrease, std::vector<double> possible_split_values, std::vector<double> &sums, std::vector<size_t> &counter);
    void findBestSplitValueLargeQ(size_t nodeID, size_t varID, double sum_node, size_t num_samples_node, double &best_value, size_t &best_varID, double &best_decrease);
    void findBestSplitValueUnordered(size_t nodeID, size_t varID, double sum_node, size_t num_samples_node, double &best_value, size_t &best_varID, double &best_decrease);

    bool findBestSplitMaxstat(size_t nodeID, std::vector<size_t> &possible_split_varIDs);

    bool findBestSplitExtraTrees(size_t nodeID, std::vector<size_t> &possible_split_varIDs);
    void findBestSplitValueExtraTrees(size_t nodeID, size_t varID, double sum_node, size_t num_samples_node, double &best_value, size_t &best_varID, double &best_decrease);
    void findBestSplitValueExtraTrees(size_t nodeID, size_t varID, double sum_node, size_t num_samples_node, double &best_value, size_t &best_varID, double &best_decrease, std::vector<double> possible_split_values, std::vector<double> &sums_right, std::vector<size_t> &n_right);
    void findBestSplitValueExtraTreesUnordered(size_t nodeID, size_t varID, double sum_node, size_t num_samples_node, double &best_value, size_t &best_varID, double &best_decrease);

    bool findBestSplitBeta(size_t nodeID, std::vector<size_t> &possible_split_varIDs);
    void findBestSplitValueBeta(size_t nodeID, size_t varID, double sum_node, size_t num_samples_node, double &best_value, size_t &best_varID, double &best_decrease);
    void findBestSplitValueBeta(size_t nodeID, size_t varID, double sum_node, size_t num_samples_node, double &best_value, size_t &best_varID, double &best_decrease, std::vector<double> possible_split_values, std::vector<double> &sums_right, std::vector<size_t> &n_right);

    void addImpurityImportance(size_t nodeID, size_t varID, double decrease);

    double computePredictionMSE();

    void cleanUpInternal() override
    {
      counter.clear();
      counter.shrink_to_fit();
      sums.clear();
      sums.shrink_to_fit();
    }

    std::vector<size_t> counter; // counter
    std::vector<double> sums; // sums 
  };

} // namespace spruce

#endif /* TREEREGRESSION_H_ */
