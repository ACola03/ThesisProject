### **Meeting Notes**

------------------------------------------------------------------------

#### **June 5, 2026:**

**Improving the Wald:**

-   Are there correction methods for Wald GLM?

-   Are there better approximations than the Wald?

-   How well do these tests work? are there rules of thumb when to use the Wald?

    -   permutation wrapper

**Plot Specifics:**

-   Statistics for piano plots

    -   Variance of bar heights

-   Getting into SlugPlots (double-check):

    -   Good visual tool for efficiency

    -   Rules of thumb about how un-normal their data should be

**Supervised and Unsupervised:**

- The original “fuzzed” piano plots are “supervised”

	- Making use of knowledge about the values of the statistic

- The piano plots we're trying with Wald are unsupervised (I think)

	- That is, based on aggregating P values, not based on knowledge about possible values

	- Try to confirm you can get perfect unsupervised piano plots for
	binomial or Poisson

-   Unsupervised fuzzing of something that is clean (exact poisson test)

    -   Then empirically fuzz for Wald when $\lambda = 50$

**Weird Wald P-Values:** (ideas that were brought up)

-   Increase $\lambda \rightarrow 50$

    -   Is 0 (count) logically 0 or 1 (p-value) ... do it logically without guessing to make my brain hurt

**Fitting Models:**

-   Fit simple regression models (1 or 2 samples)

    -   See how bad they are if assumptions are violated (residuals, non-normal response)

-   Linear Model:

    -   Estimate

    -   Not sure of assumptions

    -   Run permutations to fit model

    -   Get confidence interval from permutation $\rightarrow$ should this always be done?

**Applications:**

-   When and where do we actually use a GLM

------------------------------------------------------------------------
