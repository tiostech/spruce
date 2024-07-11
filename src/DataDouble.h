/*-------------------------------------------------------------------------------
 This file is part of spruce.

 Copyright (c) [2014-2018] [Marvin N. Wright]

 This software may be modified and distributed under the terms of the MIT license.

 Please note that the C++ core of spruce is distributed under MIT license and the
 R package "spruce" under GPL3 license.
 #-------------------------------------------------------------------------------*/

// Ignore in coverage report (not used in R package)
// #nocov start
#ifndef DATADOUBLE_H_
#define DATADOUBLE_H_

#include <vector>
#include <utility>

#include "globals.h"
#include "utility.h"
#include "Data.h"

namespace spruce {

class DataDouble: public Data {
public:
  DataDouble() = default;
  
  DataDouble(const DataDouble&) = delete;
  DataDouble& operator=(const DataDouble&) = delete;

  virtual ~DataDouble() override = default;

  double get_x(size_t row, size_t col) const override {
    // Use permuted data for corrected impurity importance
    size_t col_permuted = col;
    if (col >= num_cols) {
      col = getUnpermutedVarID(col);
      row = getPermutedSampleID(row);
    }

    if (col < num_cols_no_snp) {
      return x[col * num_rows + row];
    } else {
      return getSnp(row, col, col_permuted);
    }
  }

  double get_y(size_t row, size_t col) const override {
    return y[col * num_rows + row];
  }
  
  double get_z(size_t row, size_t col) const override{
    return z[col * num_rows + row];
  }
  
  double get_test_x(size_t row, size_t col) const override
  {
    return test_x[col * num_rows + row];
  }

  void reserveMemory(size_t y_cols) override {
    x.resize(num_cols * num_rows);
    y.resize(y_cols * num_rows);
    z.resize(y_cols * num_rows);
  }

  void set_x(size_t col, size_t row, double value, bool& error) override {
    x[col * num_rows + row] = value;
  }

  void set_y(size_t col, size_t row, double value, bool& error) override {
    y[col * num_rows + row] = value;
  }
  
  void set_z(size_t col, size_t row, double value, bool& error) override {
    z[col * num_rows + row] = value;
  }

private:
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<double> test_x;
};

} // namespace spruce

#endif /* DATADOUBLE_H_ */
// #nocov end
