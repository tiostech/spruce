/*-------------------------------------------------------------------------------
 This file is part of spruce.

 Copyright (c) [2014-2018] [Marvin N. Wright]

 This software may be modified and distributed under the terms of the MIT license.

 Please note that the C++ core of spruce is distributed under MIT license and the
 R package "spruce" under GPL3 license.
 #-------------------------------------------------------------------------------*/

#include <iterator>

#include "Tree.h"
#include "utility.h"

namespace spruce
{

  // Tree::Tree() indicates a constructor of the Tree class
  // Each line after the colon (:) initializes a specific member variable of the Tree class.
  // The values being assigned (e.g., 0, false, true, or constants like DEFAULT_SPLITRULE, DEFAULT_ALPHA, etc.) are default initial values for these member variables.
  // you cannot directly instantiate it because abstract classes are meant to be inherited by other classes, which then provide implementations for the pure virtual functions declared in the abstract class.
  Tree::Tree() : mtry(0), num_samples(0), num_samples_oob(0), min_node_size(0), min_bucket(0), deterministic_varIDs(0), split_select_weights(0), case_weights(
                                                                                                                                                     0),
                 manual_inbag(0), oob_sampleIDs(0), save_node_stats(false), num_samples_nodes(0), node_predictions(0),
                 holdout(false), keep_inbag(false), data(0), regularization_factor(0), regularization_usedepth(false),
                 split_varIDs_used(0), variable_importance(0), importance_mode(DEFAULT_IMPORTANCE_MODE), sample_with_replacement(
                                                                                                             true),
                 sample_fraction(0), memory_saving_splitting(false), splitrule(DEFAULT_SPLITRULE), alpha(DEFAULT_ALPHA), minprop(
                                                                                                                             DEFAULT_MINPROP),
                 num_random_splits(DEFAULT_NUM_RANDOM_SPLITS), max_depth(DEFAULT_MAXDEPTH), depth(0), last_left_nodeID(
                                                                                                          0)
  {
  }

  // This constructor definition for the Tree class is quite similar to the previous one but includes additional parameters that are passed as arguments to initialize the object.
  Tree::Tree(std::vector<std::vector<size_t>> &child_nodeIDs, std::vector<size_t> &split_varIDs,
             std::vector<double> &split_values) : mtry(0), num_samples(0), num_samples_oob(0), min_node_size(0), min_bucket(0), deterministic_varIDs(0), split_select_weights(0), case_weights(0), manual_inbag(0), split_varIDs(split_varIDs), split_values(split_values), child_nodeIDs(child_nodeIDs), oob_sampleIDs(0), save_node_stats(false), num_samples_nodes(0), node_predictions(0),
                                                  holdout(false), keep_inbag(false), data(0), regularization_factor(0), regularization_usedepth(false), split_varIDs_used(
                                                                                                                                                            0),
                                                  variable_importance(0), importance_mode(DEFAULT_IMPORTANCE_MODE), sample_with_replacement(true), sample_fraction(
                                                                                                                                                       0),
                                                  memory_saving_splitting(false), splitrule(DEFAULT_SPLITRULE), alpha(DEFAULT_ALPHA), minprop(
                                                                                                                                          DEFAULT_MINPROP),
                                                  num_random_splits(DEFAULT_NUM_RANDOM_SPLITS), max_depth(DEFAULT_MAXDEPTH), depth(0), last_left_nodeID(
                                                                                                                                           0)
  {
  }

  // This method is likely used to configure the Tree object before it starts building or processing data.
  // It initializes various settings and configurations that affect how the tree will be constructed or used, including parameters for splitting, sampling, and regularization.
  // Flexibility: Pointers are more flexible as they can be reassigned and can be null, whereas references cannot be re-assigned after initialization and must always refer to an object.
  // Safety: References are generally safer than pointers because they cannot be null and are always guaranteed to point to a valid object after initialization.
  // Syntax: Syntax for dereferencing and accessing members differs (* for pointers, direct access for references).
  // Typical Use:
  // Use references when you have a clear ownership or association of an object (e.g., in function parameters for passing and modifying objects).
  // Use pointers when you need more flexibility (e.g., handling arrays, dynamic memory allocation) or when dealing with potentially null objects.
  // Reference (&): Using & allows you to avoid copying large objects like vectors, which can be more efficient in terms of memory and performance.
  // In summary, both references and pointers serve similar purposes of indirect access but differ in syntax, initialization, flexibility, and safety considerations. The choice between them depends on the specific requirements of your program and the problem you are solving.

  void Tree::init(const Data *data, uint mtry, size_t num_samples, uint seed, std::vector<size_t> *deterministic_varIDs,
                  std::vector<double> *split_select_weights, ImportanceMode importance_mode, uint min_node_size, uint min_bucket,
                  bool sample_with_replacement, bool memory_saving_splitting, SplitRule splitrule, std::vector<double> *case_weights,
                  std::vector<size_t> *manual_inbag, bool keep_inbag, std::vector<double> *sample_fraction, double alpha,
                  double minprop, bool holdout, uint num_random_splits, uint max_depth, std::vector<double> *regularization_factor,
                  bool regularization_usedepth, std::vector<bool> *split_varIDs_used, bool save_node_stats)
  {

    this->data = data;
    this->mtry = mtry;
    this->num_samples = num_samples;
    this->memory_saving_splitting = memory_saving_splitting;
    this->save_node_stats = save_node_stats;

    // Create root node, assign bootstrap sample and oob samples
    child_nodeIDs.push_back(std::vector<size_t>()); // vector of vectors
    child_nodeIDs.push_back(std::vector<size_t>()); // vector of vectors
    createEmptyNode();

    // Initialize random number generator and set seed
    random_number_generator.seed(seed);

    this->deterministic_varIDs = deterministic_varIDs;
    this->split_select_weights = split_select_weights;
    this->importance_mode = importance_mode;
    this->min_node_size = min_node_size;
    this->min_bucket = min_bucket;
    this->sample_with_replacement = sample_with_replacement;
    this->splitrule = splitrule;
    this->case_weights = case_weights;
    this->manual_inbag = manual_inbag;
    this->keep_inbag = keep_inbag;
    this->sample_fraction = sample_fraction;
    this->holdout = holdout;
    this->alpha = alpha;
    this->minprop = minprop;
    this->num_random_splits = num_random_splits;
    this->max_depth = max_depth;
    this->regularization_factor = regularization_factor;
    this->regularization_usedepth = regularization_usedepth;
    this->split_varIDs_used = split_varIDs_used;

    // Regularization
    if (regularization_factor->size() > 0)
    {
      regularization = true;
    }
    else
    {
      regularization = false;
    }
  }

  // void Tree::grow(std::vector<double> *variable_importance)
  // {
  //   // Allocate memory for tree growing
  //   allocateMemory();
  // 
  //   this->variable_importance = variable_importance; // Variable importance for all variables, update after splitting
  // 
  //   // Bootstrap, dependent if weighted or not and with or without replacement
  //   if (!case_weights->empty())
  //   {
  //     if (sample_with_replacement)
  //     {
  //       bootstrapWeighted();
  //     }
  //     else
  //     {
  //       bootstrapWithoutReplacementWeighted();
  //     }
  //   }
  //   else if (sample_fraction->size() > 1)
  //   {
  //     if (sample_with_replacement)
  //     {
  //       bootstrapClassWise();
  //     }
  //     else
  //     {
  //       bootstrapWithoutReplacementClassWise();
  //     }
  //   }
  //   else if (!manual_inbag->empty())
  //   {
  //     setManualInbag();
  //   }
  //   else
  //   {
  //     if (sample_with_replacement)
  //     {
  //       bootstrap();
  //     }
  //     else
  //     {
  //       bootstrapWithoutReplacement();
  //     }
  //   }
  // 
  //   // Init start and end positions
  //   start_pos[0] = 0;              // node_ID
  //   end_pos[0] = sampleIDs.size(); // node_ID // All sampleIDs in the tree, will be re-ordered while splitting
  // 
  //   // While not all nodes terminal, split next node
  //   size_t num_open_nodes = 1;
  //   size_t i = 0;
  //   depth = 0;
  //   while (num_open_nodes > 0)
  //   {
  //     // Split node
  //     bool is_terminal_node = splitNode(i); // node_ID, split node or not 
  //     if (is_terminal_node)
  //     {
  //       --num_open_nodes;
  //     }
  //     else
  //     {
  //       ++num_open_nodes;
  //       if (i >= last_left_nodeID)
  //       {
  //         // If new level, increase depth
  //         // (left_node saves left-most node in current level, new level reached if that node is splitted)
  //         last_left_nodeID = split_varIDs.size() - 2;
  //         ++depth;
  //       }
  //     }
  //     ++i;
  //   }
  //   
  //   std::cout << ".......done grow\n";
  //   // Delete sampleID vector to save memory
  //   sampleIDs.clear();
  //   sampleIDs.shrink_to_fit();
  //   cleanUpInternal();
  // }

void Tree::grow(std::vector<double> *variable_importance)
{
  // Allocate memory for tree growing
  allocateMemory();
  
  this->variable_importance = variable_importance; // Variable importance for all variables, update after splitting
  
  // Bootstrap, dependent if weighted or not and with or without replacement
  if (!case_weights->empty())
  {
    if (sample_with_replacement)
    {
      bootstrapWeighted();
    }
    else
    {
      bootstrapWithoutReplacementWeighted();
    }
  }
  else if (sample_fraction->size() > 1)
  {
    if (sample_with_replacement)
    {
      bootstrapClassWise();
    }
    else
    {
      bootstrapWithoutReplacementClassWise();
    }
  }
  else if (!manual_inbag->empty())
  {
    setManualInbag();
  }
  else
  {
    if (sample_with_replacement)
    {
      bootstrap();
    }
    else
    {
      bootstrapWithoutReplacement();
    }
  }
  
  // Init start and end positions
  start_pos[0] = 0;              // node_ID
  end_pos[0] = sampleIDs.size(); // node_ID // All sampleIDs in the tree, will be re-ordered while splitting
  
  // While not all nodes terminal, split next node
  size_t num_open_nodes = 1;
  size_t i = 0;
  depth = 0;
  
  while (num_open_nodes > 0)
  {
    
    // find the node where the test data is in 
    // Split node
    bool is_terminal_node = splitNode(i); // node_ID, split node or not 
    if (is_terminal_node){
      --num_open_nodes;
    }
    else{
      size_t split_varID = split_varIDs[i];
      double split_value = split_values[i];
      double test_value = data -> get_test_x(0, split_varID);
      // std::cout << "split_varID = " << split_varID << ": split_value = " << split_value << ", test_value = " << test_value << "\n";
      if(test_value <= split_value){
        i = child_nodeIDs[0][i];
      }else{
        i = child_nodeIDs[1][i];
      }
      ++depth;
    }
  }
  std::cout << "depth = " << depth << ".......done growing\n";
  // Delete sampleID vector to save memory
  sampleIDs.clear();
  sampleIDs.shrink_to_fit();
  cleanUpInternal();
}


  void Tree::predict(const Data *prediction_data, bool oob_prediction)
  {
    // find terminal node
    size_t num_samples_predict;
    if (oob_prediction)
    {
      num_samples_predict = num_samples_oob;
    }
    else
    {
      num_samples_predict = prediction_data->getNumRows(); // in-sample predictions
    }

    prediction_terminal_nodeIDs.resize(num_samples_predict, 0);

    // For each sample start in root, drop down the tree and return final value
    for (size_t i = 0; i < num_samples_predict; ++i)
    {
      size_t sample_idx;
      if (oob_prediction)
      {
        sample_idx = oob_sampleIDs[i];
      }
      else
      {
        sample_idx = i;
      }
      size_t nodeID = 0;
      while (1)
      {

        // Break if terminal node
        if (child_nodeIDs[0][nodeID] == 0 && child_nodeIDs[1][nodeID] == 0)
        {
          break;
        }

        // Move to child
        size_t split_varID = split_varIDs[nodeID];

        double value = prediction_data->get_x(sample_idx, split_varID);

        if (prediction_data->isOrderedVariable(split_varID))
        {
          if (value <= split_values[nodeID])
          {
            // Move to left child
            nodeID = child_nodeIDs[0][nodeID];
          }
          else
          {
            // Move to right child
            nodeID = child_nodeIDs[1][nodeID];
          }
        }
        else
        {
          size_t factorID = floor(value) - 1;
          size_t splitID = floor(split_values[nodeID]);

          // Left if 0 found at position factorID
          if (!(splitID & (1ULL << factorID)))
          {
            // Move to left child
            nodeID = child_nodeIDs[0][nodeID];
          }
          else
          {
            // Move to right child
            nodeID = child_nodeIDs[1][nodeID];
          }
        }
      }

      prediction_terminal_nodeIDs[i] = nodeID;
      // std::cout << "nodeID = " << nodeID << "\n";
    }
  }

  void Tree::computePermutationImportance(std::vector<double> &forest_importance, std::vector<double> &forest_variance,
                                          std::vector<double> &forest_importance_casewise)
  {

    size_t num_independent_variables = data->getNumCols();

    // Compute normal prediction accuracy for each tree. Predictions already computed..
    double accuracy_normal;
    std::vector<double> prederr_normal_casewise;
    std::vector<double> prederr_shuf_casewise;
    if (importance_mode == IMP_PERM_CASEWISE)
    {
      prederr_normal_casewise.resize(num_samples_oob, 0);
      prederr_shuf_casewise.resize(num_samples_oob, 0);
      accuracy_normal = computePredictionAccuracyInternal(&prederr_normal_casewise);
    }
    else
    {
      accuracy_normal = computePredictionAccuracyInternal(NULL);
    }

    prediction_terminal_nodeIDs.clear();
    prediction_terminal_nodeIDs.resize(num_samples_oob, 0);

    // Reserve space for permutations, initialize with oob_sampleIDs
    std::vector<size_t> permutations(oob_sampleIDs);

    // Randomly permute for all independent variables
    for (size_t i = 0; i < num_independent_variables; ++i)
    {

      // Check whether the i-th variable is used in the
      // tree:
      bool isused = false;
      for (size_t j = 0; j < split_varIDs.size(); ++j)
      {
        if (split_varIDs[j] == i)
        {
          isused = true;
          break;
        }
      }

      // Only do permutations if the variable is used in the tree, otherwise variable importance is 0
      if (isused)
      {
        // Permute and compute prediction accuracy again for this permutation and save difference
        permuteAndPredictOobSamples(i, permutations);
        double accuracy_permuted;
        if (importance_mode == IMP_PERM_CASEWISE)
        {
          accuracy_permuted = computePredictionAccuracyInternal(&prederr_shuf_casewise);
          for (size_t j = 0; j < num_samples_oob; ++j)
          {
            size_t pos = i * num_samples + oob_sampleIDs[j];
            forest_importance_casewise[pos] += prederr_shuf_casewise[j] - prederr_normal_casewise[j];
          }
        }
        else
        {
          accuracy_permuted = computePredictionAccuracyInternal(NULL);
        }

        double accuracy_difference = accuracy_normal - accuracy_permuted;
        forest_importance[i] += accuracy_difference;

        // Compute variance
        if (importance_mode == IMP_PERM_BREIMAN)
        {
          forest_variance[i] += accuracy_difference * accuracy_difference;
        }
        else if (importance_mode == IMP_PERM_LIAW)
        {
          forest_variance[i] += accuracy_difference * accuracy_difference * num_samples_oob;
        }
      }
    }
  }

  // #nocov start
  void Tree::appendToFile(std::ofstream &file)
  {

    // Save general fields
    saveVector2D(child_nodeIDs, file);
    saveVector1D(split_varIDs, file);
    saveVector1D(split_values, file);

    // Call special functions for subclasses to save special fields.
    appendToFileInternal(file);
  }
  // #nocov end

  void Tree::createPossibleSplitVarSubset(std::vector<size_t> &result)
  {
    // result possible_split_VarIDs
    size_t num_vars = data->getNumCols();

    // For corrected Gini importance add dummy variables
    if (importance_mode == IMP_GINI_CORRECTED)
    {
      num_vars += data->getNumCols();
    }

    // Randomly add non-deterministic variables (according to weights if needed)
    // result has been modified
    if (split_select_weights->empty())
    {
      if (deterministic_varIDs->empty())
      {
        drawWithoutReplacement(result, random_number_generator, num_vars, mtry);
      }
      else
      {
        drawWithoutReplacementSkip(result, random_number_generator, num_vars, (*deterministic_varIDs), mtry);
      }
    }
    else
    {
      drawWithoutReplacementWeighted(result, random_number_generator, num_vars, mtry, *split_select_weights);
    }

    // Always use deterministic variables
    std::copy(deterministic_varIDs->begin(), deterministic_varIDs->end(), std::inserter(result, result.end()));
  }

  bool Tree::splitNode(size_t nodeID)
  {

    // Select random subset of variables to possibly split at
    std::vector<size_t> possible_split_varIDs; // initialize
    createPossibleSplitVarSubset(possible_split_varIDs);

    // Call subclass method, sets split_varIDs and split_values
    // bool stop = splitNodeInternalCrystal(nodeID, possible_split_varIDs);  // crystal loss 
    bool stop; 
    if(splitrule == VARIANCE){
      stop = splitNodeInternal(nodeID, possible_split_varIDs); //  split node
    }else{
      stop = splitNodeInternalV2(nodeID, possible_split_varIDs); // absolute loss
    }
    
    // return the best split_varID and split_value
    if (stop)
    {
      // Terminal node
      return true;
    }
    // not terminal node 
    // get split_varID and split_value
    size_t split_varID = split_varIDs[nodeID];
    double split_value = split_values[nodeID];

    // Save non-permuted variable for prediction
    split_varIDs[nodeID] = data->getUnpermutedVarID(split_varID);
    
    // Create child nodes
    size_t left_child_nodeID = split_varIDs.size();
    child_nodeIDs[0][nodeID] = left_child_nodeID;
    createEmptyNode();
    start_pos[left_child_nodeID] = start_pos[nodeID];

    size_t right_child_nodeID = split_varIDs.size();
    child_nodeIDs[1][nodeID] = right_child_nodeID;
    createEmptyNode();
    start_pos[right_child_nodeID] = end_pos[nodeID];
    
    // std::cout << "split_varID = " << split_varID << "\n"; 
    // std::cout << "split_value = " << split_value << "\n"; 
    // std::cout << "child_nodeIDs[0][nodeID] = " << child_nodeIDs[0][nodeID] << "\n"; 
    // std::cout << "child_nodeIDs[1][nodeID] = " << child_nodeIDs[1][nodeID] << "\n"; 
    // std::cout << "before asignment \n"; 
    // for (size_t s = 0; s < sampleIDs.size(); ++s) {
    //   std::cout << sampleIDs[s] << " ";
    // }

    // For each sample in node, assign to left or right child
    if (data->isOrderedVariable(split_varID))
    {
      // Ordered: left is <= splitval and right is > splitval
      size_t pos = start_pos[nodeID];
      while (pos < start_pos[right_child_nodeID])
      {
        size_t sampleID = sampleIDs[pos];
        if (data->get_x(sampleID, split_varID) <= split_value)
        {
          // If going to left, do nothing
          ++pos;
        }
        else
        {
          // If going to right, move to right end
          --start_pos[right_child_nodeID];
          std::swap(sampleIDs[pos], sampleIDs[start_pos[right_child_nodeID]]);
        }
      }
    }
    else
    {
      // Unordered: If bit at position is 1 -> right, 0 -> left
      size_t pos = start_pos[nodeID];
      while (pos < start_pos[right_child_nodeID])
      {
        size_t sampleID = sampleIDs[pos];
        double level = data->get_x(sampleID, split_varID);
        size_t factorID = floor(level) - 1;
        size_t splitID = floor(split_value);

        // Left if 0 found at position factorID
        if (!(splitID & (1ULL << factorID)))
        {
          // If going to left, do nothing
          ++pos;
        }
        else
        {
          // If going to right, move to right end
          --start_pos[right_child_nodeID];
          std::swap(sampleIDs[pos], sampleIDs[start_pos[right_child_nodeID]]);
        }
      }
    }
    
    // End position of left child is start position of right child
    end_pos[left_child_nodeID] = start_pos[right_child_nodeID];
    end_pos[right_child_nodeID] = end_pos[nodeID];
    
    // std::cout << "after asignment \n";
    // for (size_t s = 0; s < sampleIDs.size(); ++s) {
    //   std::cout << sampleIDs[s] << " ";
    // }
    // std::cout << "\n";
    // std::cout << "start_pos[left_child_nodeID] = " << start_pos[left_child_nodeID] << " end_pos[left_child_nodeID] = " << end_pos[left_child_nodeID] << "\n";
    // std::cout << "start_pos[right_child_nodeID] = " << start_pos[right_child_nodeID] << " end_pos[right_child_nodeID] = " << end_pos[right_child_nodeID] << "\n";
    // 
    
    // No terminal node
    return false;
  }

  void Tree::createEmptyNode()
  {
    split_varIDs.push_back(0);     // node_IDs
    split_values.push_back(0);     // node_IDs
    child_nodeIDs[0].push_back(0); // node_ID
    child_nodeIDs[1].push_back(0); // node_ID
    start_pos.push_back(0);        // node_ID
    end_pos.push_back(0);          // node_ID

    if (save_node_stats)
    {
      num_samples_nodes.push_back(0); // node_ID
      split_stats.push_back(0);       // node_ID
    }

    createEmptyNodeInternal(); // node_ID
  }

  size_t Tree::dropDownSamplePermuted(size_t permuted_varID, size_t sampleID, size_t permuted_sampleID)
  {

    // Start in root and drop down
    size_t nodeID = 0;
    while (child_nodeIDs[0][nodeID] != 0 || child_nodeIDs[1][nodeID] != 0)
    {

      // Permute if variable is permutation variable
      size_t split_varID = split_varIDs[nodeID];
      size_t sampleID_final = sampleID;
      if (split_varID == permuted_varID)
      {
        sampleID_final = permuted_sampleID; // permuted....
      }

      // Move to child
      double value = data->get_x(sampleID_final, split_varID);
      if (data->isOrderedVariable(split_varID))
      {
        if (value <= split_values[nodeID])
        {
          // Move to left child
          nodeID = child_nodeIDs[0][nodeID];
        }
        else
        {
          // Move to right child
          nodeID = child_nodeIDs[1][nodeID];
        }
      }
      else
      {
        size_t factorID = floor(value) - 1;
        size_t splitID = floor(split_values[nodeID]);

        // Left if 0 found at position factorID
        if (!(splitID & (1ULL << factorID)))
        {
          // Move to left child
          nodeID = child_nodeIDs[0][nodeID];
        }
        else
        {
          // Move to right child
          nodeID = child_nodeIDs[1][nodeID];
        }
      }
    }
    return nodeID;
  }

  void Tree::permuteAndPredictOobSamples(size_t permuted_varID, std::vector<size_t> &permutations)
  {

    // Permute OOB sample
    // std::vector<size_t> permutations(oob_sampleIDs);
    std::shuffle(permutations.begin(), permutations.end(), random_number_generator);

    // For each sample, drop down the tree and add prediction
    // outofbag permuated
    for (size_t i = 0; i < num_samples_oob; ++i)
    {
      size_t nodeID = dropDownSamplePermuted(permuted_varID, oob_sampleIDs[i], permutations[i]);
      prediction_terminal_nodeIDs[i] = nodeID;
    }
  }

  void Tree::bootstrap()
  {

    // Use fraction (default 63.21%) of the samples
    size_t num_samples_inbag = (size_t)num_samples * (*sample_fraction)[0];

    // Reserve space, reserve a little more to be save)
    sampleIDs.reserve(num_samples_inbag);
    oob_sampleIDs.reserve(num_samples * (exp(-(*sample_fraction)[0]) + 0.1));

    std::uniform_int_distribution<size_t> unif_dist(0, num_samples - 1);

    // Start with all samples OOB
    inbag_counts.resize(num_samples, 0);

    // Draw num_samples samples with replacement (num_samples_inbag out of n) as inbag and mark as not OOB
    for (size_t s = 0; s < num_samples_inbag; ++s)
    {
      size_t draw = unif_dist(random_number_generator);
      sampleIDs.push_back(draw);
      ++inbag_counts[draw];
    }

    // Save OOB samples
    for (size_t s = 0; s < inbag_counts.size(); ++s)
    {
      if (inbag_counts[s] == 0)
      {
        oob_sampleIDs.push_back(s);
      }
    }
    num_samples_oob = oob_sampleIDs.size();

    if (!keep_inbag)
    {
      inbag_counts.clear();
      inbag_counts.shrink_to_fit();
    }
  }

  void Tree::bootstrapWeighted()
  {

    // Use fraction (default 63.21%) of the samples
    size_t num_samples_inbag = (size_t)num_samples * (*sample_fraction)[0];

    // Reserve space, reserve a little more to be save)
    sampleIDs.reserve(num_samples_inbag);
    oob_sampleIDs.reserve(num_samples * (exp(-(*sample_fraction)[0]) + 0.1));

    std::discrete_distribution<> weighted_dist(case_weights->begin(), case_weights->end());

    // Start with all samples OOB
    inbag_counts.resize(num_samples, 0);

    // Draw num_samples samples with replacement (n out of n) as inbag and mark as not OOB
    for (size_t s = 0; s < num_samples_inbag; ++s)
    {
      size_t draw = weighted_dist(random_number_generator);
      sampleIDs.push_back(draw);
      ++inbag_counts[draw];
    }

    // Save OOB samples. In holdout mode these are the cases with 0 weight.
    if (holdout)
    {
      for (size_t s = 0; s < (*case_weights).size(); ++s)
      {
        if ((*case_weights)[s] == 0)
        {
          oob_sampleIDs.push_back(s);
        }
      }
    }
    else
    {
      for (size_t s = 0; s < inbag_counts.size(); ++s)
      {
        if (inbag_counts[s] == 0)
        {
          oob_sampleIDs.push_back(s);
        }
      }
    }
    num_samples_oob = oob_sampleIDs.size();

    if (!keep_inbag)
    {
      inbag_counts.clear();
      inbag_counts.shrink_to_fit();
    }
  }

  void Tree::bootstrapWithoutReplacement()
  {

    // Use fraction (default 63.21%) of the samples
    size_t num_samples_inbag = (size_t)num_samples * (*sample_fraction)[0];
    shuffleAndSplit(sampleIDs, oob_sampleIDs, num_samples, num_samples_inbag, random_number_generator);
    num_samples_oob = oob_sampleIDs.size();

    if (keep_inbag)
    {
      // All observation are 0 or 1 times inbag
      inbag_counts.resize(num_samples, 1);
      for (size_t i = 0; i < oob_sampleIDs.size(); i++)
      {
        inbag_counts[oob_sampleIDs[i]] = 0;
      }
    }
  }

  void Tree::bootstrapWithoutReplacementWeighted()
  {

    // Use fraction (default 63.21%) of the samples
    size_t num_samples_inbag = (size_t)num_samples * (*sample_fraction)[0];
    drawWithoutReplacementWeighted(sampleIDs, random_number_generator, num_samples - 1, num_samples_inbag, *case_weights);

    // All observation are 0 or 1 times inbag
    inbag_counts.resize(num_samples, 0);
    for (auto &sampleID : sampleIDs)
    {
      inbag_counts[sampleID] = 1;
    }

    // Save OOB samples. In holdout mode these are the cases with 0 weight.
    if (holdout)
    {
      for (size_t s = 0; s < (*case_weights).size(); ++s)
      {
        if ((*case_weights)[s] == 0)
        {
          oob_sampleIDs.push_back(s);
        }
      }
    }
    else
    {
      for (size_t s = 0; s < inbag_counts.size(); ++s)
      {
        if (inbag_counts[s] == 0)
        {
          oob_sampleIDs.push_back(s);
        }
      }
    }
    num_samples_oob = oob_sampleIDs.size();

    if (!keep_inbag)
    {
      inbag_counts.clear();
      inbag_counts.shrink_to_fit();
    }
  }

  void Tree::bootstrapClassWise()
  {
    // Empty on purpose (virtual function only implemented in classification and probability)
  }

  void Tree::bootstrapWithoutReplacementClassWise()
  {
    // Empty on purpose (virtual function only implemented in classification and probability)
  }

  void Tree::setManualInbag()
  {
    // Select observation as specified in manual_inbag vector
    sampleIDs.reserve(manual_inbag->size());
    inbag_counts.resize(num_samples, 0);
    for (size_t i = 0; i < manual_inbag->size(); ++i)
    {
      size_t inbag_count = (*manual_inbag)[i];
      if ((*manual_inbag)[i] > 0)
      {
        for (size_t j = 0; j < inbag_count; ++j)
        {
          sampleIDs.push_back(i);
        }
        inbag_counts[i] = inbag_count;
      }
      else
      {
        oob_sampleIDs.push_back(i);
      }
    }
    num_samples_oob = oob_sampleIDs.size();

    // Shuffle samples
    std::shuffle(sampleIDs.begin(), sampleIDs.end(), random_number_generator);

    if (!keep_inbag)
    {
      inbag_counts.clear();
      inbag_counts.shrink_to_fit();
    }
  }

} // namespace spruce
