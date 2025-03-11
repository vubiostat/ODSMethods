# ODSMethods

This package is a draft at present to bring together a wide variety of outcome dependent sampling (ODS) methods led by professors at Vanderbilt University Medical Center. 

## Healthy Package Development

* All test cases PASS on before every push to main repository.
* CRAN Check should have 0 ERRORs and 0 WARNINGS before every push to the main repository. This prevents a huge amount of work later.
* No files added to repository unless they need to be. This is intended to be a working package _not a pile._
* Work incrementally. Change one or add one small thing at a time. 
* User interface is super important at this stage, think carefully about a user who knows little about these tools and will just take data and play. Make discovery easy. Follow the principal of least surprise.
* Never use a "." dot in a function or variable name. This is a very strong recommendation from the core R team as it can cause problems with S3 dispatch--which this package will rely on.
* Reference code goes into a "helper" file in tests. This is code that has been validated via a published study. It should _**not**_ be edited in any way.
* Routines in the "R" directory are the working published and polished code.
* A routine in the "R" directory needs a tedious amount of tests that consider not only user interface, but the results and compares them with reference code.
* Proper author reference/license in each code file, and the author needs to be part of the DESCRIPTION.

## Needs

* Go broad and outline what every method would look like in the current interface approach. Right now it's limited to acml and continuous logitudinal designs. 
* Go deep on the existing acml routine and consider all the output possible using S3 (print, format, summary, plot, vcov, coef, resid). Use the existing 'lm' as a guide. 
* Start optimization of acml underlying algorithms.
  * Get rid of `nlm`. This is not a recommended optimizer.
  * Add in the later version of ACML.LME that include BLUP and other methods.
  * Add in Lucy's analytical Hessian. 
  * Last minute: Push things to C++. 
* Some reference data sets to write examples from (with full documentation). Existing reference set is really a test set and should be moved to `inst` and loaded from there. 


