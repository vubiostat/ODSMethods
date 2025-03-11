  ##############################################################################
 #
# ODSMethods Statistical methods in outcome dependent sampling
#
# Copyright (C) 2025 Shawn P. Garbett, Bailu Yan
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>

#' Fit model using ascertainment corrected likelihood model (ACML)
#'
#' Outcome dependent sampling designs need to be corrected when fitting a
#' statistical model for proper inferences. This routine will fit and return
#' the model fit using ascertainment corrected likelihood.
#'
#' @param formula `formula`; (or an object that can be
#'   coerced to that class): a symbolic description of the model to be fitted.
#'   The details of model specification are given under ‘Details’.
#' @param design an object of class 'odsdesign'. This specifies the
#'   design of sampling used for the fitting algorithm.
#' @param data `data.frame`; an optional data frame, list or environment (or
#'   object coercible by as.data.frame to a data frame) containing the variables
#'   in the model. If not found in data, the variables are taken from
#'   environment(formula), typically the environment from which acml is called.
#' @param subset an optional vector specifying a subset of observations to be
#'   used in the fitting process. (See additional details about how this
#'   argument interacts with data-dependent bases in the ‘Details’ section of
#'   the model.frame documentation.)
#' @param weights	an optional vector of weights to be used in the fitting
#'   process. Should be NULL or a numeric vector. If non-NULL, weighted least
#'   squares is used with weights weights (that is, minimizing sum(w*e^2));
#'   otherwise ordinary least squares is used. See also ‘Details’,
#' @param na.action a function which indicates what should happen when the data contain NAs. The default is set by the na.action setting of options, and is na.fail if that is unset. The ‘factory-fresh’ default is na.omit. Another possible value is NULL, no action. Value na.exclude can be useful.
#' @param subset `logical`; an optional vector specifying a subset of
#'   observations to be used in the fitting process. (See additional details
#'   about how this argument interacts with data-dependent bases in the
#'   ‘Details’ section of the model.frame documentation.)
#' @param na.action `function`; an optional vector specifying a subset of
#'   observations to be used in the fitting process. (See additional details
#'   about how this argument interacts with data-dependent bases in the
#'   ‘Details’ section of the model.frame documentation.)
#' @param verbose `numeric(1)`; Debugging information printing level
#'   (sent to the optimizer).
#' @param init `numeric`; Initial starting position for parameter search
#'   via the optimizer.
#' @param ... Optional additional parameters passed to sub methods.
#'
#' @importFrom stats lm
#' @export
acml <- function(
  formula,
  design,
  data = NULL,
  subset = NULL,
  weights = NULL,
  na.action = getOption('na.action'),
  verbose = 0L,
  init = NULL,
  ...)
{
  NULL
}
