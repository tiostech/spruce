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

##' @export
timepoints <- function(x, ...)  UseMethod("timepoints")

##' Extract unique death times of spruce Survival prediction object.
##'
##'
##' @title spruce timepoints
##' @param x spruce Survival prediction object.
##' @param ... Further arguments passed to or from other methods.
##' @return Unique death times
##' @seealso \code{\link{spruce}}
##' @author Marvin N. Wright
##' @export
timepoints.spruce.prediction <- function(x, ...) {
  if (!inherits(x, "spruce.prediction")) {
    stop("Object ist no spruce.prediction object.")
  }
  if (x$treetype != "Survival") {
    stop("No timepoints found. Object is no Survival prediction object.")
  }
  if (is.null(x$unique.death.times)) {
    stop("No timepoints found.")
  }
  return(x$unique.death.times)
}

##' Extract unique death times of spruce Survival forest
##'
##'
##' @title spruce timepoints
##' @param x spruce Survival forest object.
##' @param ... Further arguments passed to or from other methods.
##' @return Unique death times
##' @seealso \code{\link{spruce}}
##' @author Marvin N. Wright
##' @aliases timepoints
##' @export
timepoints.spruce <- function(x, ...) {
  if (!inherits(x, "spruce")) {
    stop("Object ist no spruce object.")
  }
  if (x$treetype != "Survival") {
    stop("No timepoints found. Object is no Survival forest.")
  }
  if (is.null(x$unique.death.times)) {
    stop("No timepoints found.")
  }
  return(x$unique.death.times)
}
