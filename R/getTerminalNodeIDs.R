# -------------------------------------------------------------------------------
#   This file is part of spruce.
#
# spruce is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# spruce is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with spruce. If not, see <http://www.gnu.org/licenses/>.
#
# Written by:
#
#   Marvin N. Wright
# Institut fuer Medizinische Biometrie und Statistik
# Universitaet zu Luebeck
# Ratzeburger Allee 160
# 23562 Luebeck
# Germany
#
# http://www.imbs-luebeck.de
# -------------------------------------------------------------------------------

##' This function is deprecated. 
##' Please use predict() with \code{type = "terminalNodes"} instead.
##' This function calls predict() now. 
##'
##' @title Get terminal node IDs (deprecated)
##' @param rf \code{spruce} object.
##' @param dat New dataset. Terminal node IDs for this dataset are obtained. 
##'
##' @return Matrix with terminal nodeIDs for all observations in dataset and trees.
##'
##' @examples
##' rf <- spruce(Species ~ ., data = iris, num.trees = 5, write.forest = TRUE)
##' getTerminalNodeIDs(rf, iris)
##' @export
getTerminalNodeIDs <- function(rf, dat) {
  warning("Function getTerminalNodeIDs() deprecated, calling predict().")
  predict(rf, dat, type = "terminalNodes")$predictions
}

