# ODSMethods

This package is a draft at present to bring together a wide variety of outcome dependent sampling (ODS) methods led by professors at Vanderbilt University.

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
* It should be dependency adverse. An imports for say tidyverse creates a huge dependency liability. Suggests for a package is okay, but to be avoided and should include required checks that a package is loaded if needed.

## Needs

* Need a slide show of goals / design / result.
* [Deep] Finish acml to every deep corner and usage.
  * coef() should return transformed coefficients (test should be raw). Maybe coef(model, raw=TRUE)?
  * Add names to R estimates including the 4 additional parameters
  * Add vcov that returns both levels in the manner of lme4
  * Add computation of rank/df for variables
  * Add 'fitted.values' to acml
  * Add 'residuals' to acml (and S3 routine)
  * Add format / summary to acml, coefficients, confidence intervals, p-values.
  * Add plot to acml
  * Add user interface tests to acml
  * Add sample method to odsdesign
  NOTE: In general use lm() and lmer() object output and functions as guide.  
* [Refine] acml. 
  * Replace with ACML.LME (not the validated test version!)
  * Add the analytical Hessian (Lucy)
  * Optimize for speed.
* [Document]
  * Add reference data sets preferably from papers.
  * Add vignette on usage. (start with acml)
* [Expand] Fill in wide and deep for all routines identified in specification.
* [Get Credit] Write paper for Journal of Statistical Software


### TO DO LIST

#### Lucy

* Read a few _more_ sections of "Writing R Extensions"
  ~~Lucy read 3 sections and will continue~~
* Create 3 slides (include title and credits). 1 slide with goals/example/?
* Get SMLE running, in the linear model settings
* Add format / summary to acml, "print.acml function"

#### Shawn

* Investigate convergence issues with BLAS differences. specifically vcov is way off on mac.
  Apple M1 it fails on.
  https://www.intel.com/content/www/us/en/developer/articles/technical/introduction-to-the-conditional-numerical-reproducibility-cnr.html
* Add raw=TRUE to coef, vcov

### Notes on lm S3 Methods

```
> print(sloop::s3_methods_class("lm"), n=40)
# A tibble: 40 × 4
   generic        class visible source             
   <chr>          <chr> <lgl>   <chr>              
 1 add1           lm    FALSE   registered S3method
 2 addterm        lm    FALSE   registered S3method
 3 alias          lm    FALSE   registered S3method
 4 anova          lm    FALSE   registered S3method
 5 boxcox         lm    FALSE   registered S3method
 6 case.names     lm    FALSE   registered S3method
 7 confint        lm    TRUE    stats              
 8 cooks.distance lm    FALSE   registered S3method
 9 deviance       lm    FALSE   registered S3method
10 dfbeta         lm    FALSE   registered S3method
11 dfbetas        lm    FALSE   registered S3method
12 drop1          lm    FALSE   registered S3method
13 dropterm       lm    FALSE   registered S3method
14 dummy.coef     lm    TRUE    stats              
15 effects        lm    FALSE   registered S3method
16 extractAIC     lm    FALSE   registered S3method
17 family         lm    FALSE   registered S3method
18 formula        lm    FALSE   registered S3method
19 hatvalues      lm    FALSE   registered S3method
20 influence      lm    FALSE   registered S3method
21 kappa          lm    TRUE    base               
22 labels         lm    FALSE   registered S3method
23 logLik         lm    FALSE   registered S3method
24 logtrans       lm    FALSE   registered S3method
25 model.frame    lm    FALSE   registered S3method
26 model.matrix   lm    TRUE    stats              
27 nobs           lm    FALSE   registered S3method
28 plot           lm    FALSE   registered S3method
29 predict        lm    TRUE    stats              
30 print          lm    FALSE   registered S3method
31 proj           lm    FALSE   registered S3method
32 qqnorm         lm    FALSE   registered S3method
33 qr             lm    FALSE   registered S3method
34 residuals      lm    TRUE    stats              
35 rstandard      lm    FALSE   registered S3method
36 rstudent       lm    FALSE   registered S3method
37 simulate       lm    FALSE   registered S3method
38 summary        lm    TRUE    stats              
39 variable.names lm    FALSE   registered S3method
40 vcov           lm    FALSE   registered S3method
```

it probably makes sense to add all of them, even if you have a stop message that says "this method isn't implemented yet" -- Cole Beck 3/12/25

## Example Updated Calls to Other Methods

To include BLUP, 

ods(..., BLUP="bivariate") (or "slope", or "intercept") is all that would be required. Default BLUP=NULL. 

This also requires the updated Chiara code integrated.

Thinking about SIEVE

* Y is fine, comes from _method_(formula)
* X is fine, comes from _method_(formula)
* Z would come from _design_(formula) call

What is Delta? Survival event, check survival packages how this is specified.

Bspline_Z should come from the formula of the _method_.
scale, maxiter, tol,  are control values.

His model value, "linear", "logistic" or "coxph". 

_method_ would be SMLE in this case.

