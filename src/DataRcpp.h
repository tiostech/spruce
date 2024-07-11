/*-------------------------------------------------------------------------------
 This file is part of spruce.
 spruce is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

spruce is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with spruce. If not, see <http://www.gnu.org/licenses/>.

Written by:

Marvin N. Wright
Institut für Medizinische Biometrie und Statistik
Universität zu Lübeck
Ratzeburger Allee 160
23562 Lübeck

http://www.imbs-luebeck.de
#-------------------------------------------------------------------------------*/

#ifndef DATARCPP_H_
#define DATARCPP_H_

#include <Rcpp.h>

#include "globals.h"
#include "utility.h"
#include "Data.h"

namespace spruce
{

  class DataRcpp : public Data
  {
  public:
    DataRcpp() = default;
    DataRcpp(Rcpp::NumericMatrix &x, Rcpp::NumericMatrix &y, Rcpp::NumericMatrix &z, Rcpp::NumericMatrix &test_x, std::vector<std::string> variable_names, size_t num_rows, size_t num_cols)
    {
      this->x = x;
      this->y = y;
      this->z = z;
      this->test_x = test_x;
      this->variable_names = variable_names;
      this->num_rows = num_rows;
      this->num_cols = num_cols;
      this->num_cols_no_snp = num_cols;
    }

    DataRcpp(const DataRcpp &) = delete;
    DataRcpp &operator=(const DataRcpp &) = delete;

    virtual ~DataRcpp() override = default;

    double get_x(size_t row, size_t col) const override
    {
      // Use permuted data for corrected impurity importance
      size_t col_permuted = col;
      if (col >= num_cols)
      {
        col = getUnpermutedVarID(col);
        row = getPermutedSampleID(row);
      }

      if (col < num_cols_no_snp)
      {
        return x(row, col);
      }
      else
      {
        return getSnp(row, col, col_permuted);
      }
    }

    double get_y(size_t row, size_t col) const override
    {
      return y(row, col);
    }
    double get_z(size_t row, size_t col) const override
    {
      return z(row, col);
    }
    
    double get_test_x(size_t row, size_t col) const override
    {
      return test_x(row, col);
    }

    // #nocov start
    void reserveMemory(size_t y_cols) override
    {
      // Not needed
    }

    void set_x(size_t col, size_t row, double value, bool &error) override
    {
      x(row, col) = value;
    }

    void set_y(size_t col, size_t row, double value, bool &error) override
    {
      y(row, col) = value;
    }

    void set_z(size_t col, size_t row, double value, bool &error) override
    {
      z(row, col) = value;
    }
    // #nocov end

  private:
    Rcpp::NumericMatrix x;
    Rcpp::NumericMatrix y;
    Rcpp::NumericMatrix z;
    Rcpp::NumericMatrix test_x; 
  };

} // namespace spruce

#endif /* DATARCPP_H_ */
