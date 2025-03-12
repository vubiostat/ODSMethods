# ODSMethods

This package is a draft at present to bring together a wide variety of outcome dependent sampling (ODS) methods led by professors at Vanderbilt University Medical Center. 

## Healthy Package Development

* All test cases PASS before every push to main repository.
* CRAN Check should have 0 ERRORs and 0 WARNINGS before every push to the main repository. This prevents a huge amount of work later.
* No files added to repository unless they need to be. This is intended to be a working package _not a pile._
* Work incrementally. Change one or add one small thing at a time. 
* User interface is super important at this stage, think carefully about a user who knows little about these tools and will just take data and play. Make discovery easy. Follow the principal of least surprise.
* Never use a "." dot in a function or variable name. This is a very strong recommendation from the core R team as it can cause problems with S3 dispatch--which this package will rely on.
* Reference code goes into a "helper" file in tests. This is code that has been validated via a published study. It should _**not**_ be edited in any way.
* Routines in the "R" directory is the working published and polished code.
* A routine in the "R" directory needs a tedious amount of tests that consider not only user interface, but the results and compares them with reference code.
* Proper author reference/license in each code file, and the author needs to be part of the DESCRIPTION.

## Needs

* [Wide, specification] Generate example calls for every method under consideration.
* [Deep] Finish acml to every deep corner and usage.
  * Add all references and edit R/ODSMethods-package.R
  * Add format / summary to ods
  * Add format / summary to acml
  * Add names to R estimates including the 4 additional parameters
  * Add vcov that returns both levels in the manner of lme4
  * Add computation of rank/df for variables
  * Add 'fitted.values' to acml
  * Add 'residuals' to acml (and S3 routine)
  * Add plot to acml
  NOTE: In general use lm() and lmer() object output and functions as guide.  
* [Refine] acml. 
  * Replace with ACML.LME (not the validated test version!)
  * Add the analytical Hessian.
  * Get rid of `nlm`. This is not a recommended optimizer.
  * Optimize for speed.
* [Document]
  * Add reference data sets preferably from papers.
  * Add vignette on usage. (start with acml)
* [Expand] Fill in wide and deep for all routines identified in specification.
* [Get Credit] Write paper for Journal of Statistical Software

