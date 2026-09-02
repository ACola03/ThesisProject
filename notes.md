**Meeting Notes: August 13th**

1.  coverage simulations: show coverage is between 4 and 6%; test the tail coverage separately, how big do our samples need to be before we expect to reliably have CI between 2 and 3% at each tail
2.  slugplots for the mean and shape parameter since they seem uncorrelated
3.  quantile based intervals are scale independent
    i.  HD based interval are not scale independent
    ii. frequentist example of this \^
    iii. looking at something in two directions is analogous to a HPD with alpha/2 and is scale dependent; looking at something in one direction is analogous to a quantile approach with alpha and is scale independent. If you have 1 day, asking if 2 days or 3 days is far away will always give the same answer regardless of the scale. If you look in two directions instead and ask whether 0.25days or 2 days is farther, the answer now depends on the scale of time or rate.
4.  misc
    i.  make the grey dots easier to see for JD
    ii. put histograms under the contours

Other:

1.  histogram gives more information about coverage
2.  do coverage in a simpler way using a binomial confidence interval
3.  how many times we are below 2.5% and above 97.5%

------------------------------------------------------------------------

**Meeting Notes: August 5th**

**Lab Meeting Presentation:**

1.  Better Introduction:
    i.  Describe it as a method to evaluate statistical tests.
2.  Plot Related:
    i.  Should specify that the plots on slide 5 are samples that are all generated from a $\mathrm{N}(5,1)$ distribution, but the p-values are computed by modifying the mean to 4, 5, and 6.
    ii. It might make sense to have Super Wavy appear before Wavy, since it can mistakenly be uniform from the statistics, but visually isn't.
    iii. The Violin plots should have the y-axis scaled and have one panel shown per page. Also, the violins which achieve larger statistics should be changed to have different Beta parameters so that there isn't as much of a skew.
    iv. The CvM and AD visualization should be done manually for one-sample.
3.  Code Related:
    i.  Don't hard-code the colors and line types for the power analysis plots. Even though the structure of the plots are fixed, it should be done more naturally.
4.  Misc:
    i.  Explain why the minimum is used in the calculation of a two-sided p-value. We don't care about direction, we just want to know a general understanding of extremeness.
    ii. We explored Symmetric Fuzzing but I should make the code do this form of fuzzing instead of left fuzzing. It isn't completely accurate to not "not yet fully explored".
    iii. The fourth point on slide 22 doesn't need to be there; the following point explains it better.
    iv. There should be an introduction for order statistics on slide 16. Also, the Beta argument should be written under the Space and Location bullets to be more clear.
    v.  Use "drawbacks" instead of "downfalls"

**AD Report: ✓**

1.  Show AD for different Beta
2.  Confidence value on an AD statistic
3.  Properly make the one-sample visualization
4.  Contrast it to CvM to show why tail weighting matters.

**On Fuzzing:**

1.  Could you fuzz first then flip? Instead of flipping then fuzzing then converting back.
2.  If we are validating our method, there is no reason to count ties against ourselves
3.  For empirical, we assume the methodology is conservative, we try to work against that
4.  In general, we are doing empirical fuzzing - if we do exact fuzzing, we can get exact pianos

**On Statistic Validation (the curves): ✓**

1.  CvM for non-small sample sizes
2.  CvM validity for all shapes
3.  For a fixed value of the shape ( of 3/(a+b+1) ), is it always minimized for a mean of 0.5
4.  If we have a shift, maybe we don't really care about a shape of 1 (?)

**Further Testing: Beta MLE ✓**

1.  Should we be estimating beta parameters? Do an MLE fit of a Beta to a p-value distribution. What parameters are good enough, and how many deviates would it take to say that we are confident that we are in that range.
2.  Based on the point above, at some point come up with tests that aren't good enough.

**On HPD: quantile.R**

1.  Example on poisson waiting times and concept of scale invariance. The scales we transformed from rates to waiting times (inversion). Then intervals were constructed independently on the two different scales. However, when constructing the interval on the rate then transforming the scale into waiting time, the intervals were different compared to if it had been only done on waiting time.
2.  In a poisson process you have a unit of time, and observe say 4 events (which is usually just exact (?)). Do a non-exact with MLE approach and get noticeably different CI than if you do the exact quantile approach.

**On Equivalence Tests:**

1.  This is in relation to Ben's email on 'Equivalence Tests for Piano Plots?'
2.  Since this is an equivalence test, we are caring about effect size rather than strictly p-values. For example, setting a threshold beforehand on 'what is small' and 'what is large'.
3.  If I recall correctly, it is harder to work with effect sizes instead of p-values, but I cannot think of the reasoning behind it.

------------------------------------------------------------------------

##### 2026 Jul 08 (Wed)

Adam to send link to a nicer version of provisional notes

We have many uniformity statistics: \* Histogram variance (requires boxing) \* variance of spaces (analytically nice) \* start at 0, end at 1 for n+1 spaces \* naive location \* weighted location

The location statistics seem terrible! This is a big surprise and we should probably try to understand it at least for a while. In the meantime, let's go forward with the spacing statistic.

##### 2026 Jul 03 (Fri)

Order-statistic estimators. Currently liking [the WLSE from here](https://arxiv.org/pdf/2107.09316):

------------------------------------------------------------------------

#### **Week of June 12, 2026:**

##### Meeting

We want to compare distributions of Wald P values with distributions of LRT P values (maybe obtainable through anova())

The simplest statistic to try for piano plots is the variance of the number of observations per box; is there something analogous though where we don't have to choose a box size?

We should have some idea of the null distribution of this statistic for different values of N and b (number of p values and number of boxes)

You can decide how to correct the Hauck-Donners to 0 or 1 based on the behaviour of your P value procedure in that neighborhood (i.e., the P value you are using for x=1, which is the next-door neighbor).

We have a possible new way of thinking about supervised fuzzing. In supervised fuzzing, we know the probabilities that the p-values are based on. Does that work in the case of unsupervised fuzzing?

I think it might, but we're definitely going to have to set the Hauck- Donner aside while we try to figure it out. So for now, let's try to focus on the Poisson test followed by the likelihood ratio test. We probably shouldn't rely at all on the density of observations to fuzz things, but we might be allowed to rely on the nominal p-values that we're getting. We probably have to start by assuming that the p-values are conservative, since that's what p-values are supposed to be. And that kind of breaks the idea of handling Hauck-Donner as a borderline value.

Start by reading what Roswell wrote about the fuzzing logic. See if you think it's clear and either quote it where we can find it or see if you can make it more clear.

------------------------------------------------------------------------

#### **Week of June 5, 2026:**

##### Meeting

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

-   The original “fuzzed” piano plots are “supervised”

-   Making use of knowledge about the values of the statistic

-   The piano plots we're trying with Wald are unsupervised (I think)

-   That is, based on aggregating P values, not based on knowledge about possible data values

-   Try to confirm you can get perfect unsupervised piano plots for binomial or Poisson

-   Unsupervised fuzzing of something that is clean (exact poisson test)

-   Then empirically fuzz for Wald when $\lambda = 50$

-   and then for LRT

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

\newpage

#### **June 9, 2026 (pre-meeting thoughts about the above):**

**Improving The Wald:**

-   Correction methods:

-   I haven't come up with anything yet, but the Wald normal approximation assigns a lower p-value to the mean response than the exact test.

-   For example, if $\lambda=3$, the exact p-value is $0.6472$ whereas the Wald gives $0.5$. So this might be the reason why the piano plots I was showing last week were slightly left-biased towards smaller p-values. Knowing this, a correction method might involve pushing these in the opposite direction ... but then it's no longer a Wald test (which should be perfectly fine).

-   Poisson calculator link when I don't want to use R: <https://stattrek.com/online-calculator/poisson>

-   Better approximations than the Wald:

-   It appears that Likelihood-Ratio Tests and Score Tests are existing alternatives

-   Links (*better approximations than the wald for poisson family glm*): [link1](https://stats.stackexchange.com/questions/449344/likelihood-ratio-wald-and-score-are-equivalent){.uri}

-   Google AI refers to the issue as the **Hauck–Donner effect**, where the standard error inflates faster than the coefficient, leading to falsely non-significant results for large coefficients.

-   I will look into these once we understand *Supervised* and *Unsupervised*

**Plot Specifics:**

-   Statistics for pianoPlots:

-   When looking at the variance of bar heights, we *obviously* prefer a small variance to indicate that the p-values are uniformly spread across the domain. However, should we come up with rules of thumb on how variable the heights should be? or leave it as a clarity interpretation where we don't specifically set any threshold and otherrwise let the user decide whether the variance is large enough to make it unclear whether the p-values uniformly distributed.

-   Getting into slugPlots:

-   The slugPlots are where we order the confidence intervals and ideally, 2.5% are too low, and 2.5% are too high (missing the true value of the estimate).

-   I still have to look into how the confidence intervals are ordered, but an introduction is given in *Assessing Hill Diversity.*

-   When we said *good tool for efficiency*, was the implied meaning that it allows us to see how often our confidence intervals are too low/high (the purpose of the plot)

**Supervised and Unsupervised:**

-   The original “fuzzed” piano plots are “supervised”:

-   If we say this, then are we implying that using p-values from the exact null distribution (not any approximations) is a supervised approach.

-   If we assume the above, I cannot then apply an unsupervised fuzzing approach to something that is clean.

-   So, I think it might make more sense to have the distinction of supervision defined as making use of the poisson count that generated that p-value. To be specific, in the Wald case, the count of 0 is binned together with the mean response, which shouldn't occur if 0 has some probability mass.

-   If we are unsupervised, we only look at the obtained wald p-values, ignore the counts that generated them, order the p-values, and fuzz within the lagged intervals. If we are supervised, we fuzz based on the ordering of the count instead of the p-value; an issue I found is that due to the inflated p-value at 0, some intervals will have a minimum larger than the maximum ... (see pdf from last meeting)

-   Is this a correct interpretation of supervision?

**Weird Wald P-Values:**

-   Increase $\lambda \rightarrow50$ results sent by email (since they won't show on github repo):
-   Non-Fuzzed Wald P-Values
-   Fuzzed P-Values - excluding any counts of 0 (1,10,50)
-   Fuzzed P-Values - including all counts of 0 (1,10,50) ... main object of interest
-   As expected, larger $\lambda$ produces better pianos

**Fitting Models:**

-   These points were not yet explored since the above is a higher priority.

**Applications:**

-   Not explored, but can be discussed if time permits.

##### 2026 Jul 22 (Wed)

Results are interesting but also kind of disappointing. We don't have any test that produces a very tight distribution with 1,000 samples. We had decided to check 2,000 samples and Adam should do that, but it's not likely to help much.

DTS looks like it's a two sample test. Adam is going to look into whether there's a one-sample state-of-the-art location-based test and if it's equivalent to Anderson-Darling, and implement it.

Give some thought to how you normalize things. We should understand better the behavior of the statistics with respect to little n, the size of the sample. If things are normalized, the mean shouldn't be changing and we should be able to tell if the distribution is getting narrower or not. We had some questions about this bill. Excuse me. We had some questions about this with the location statistic before.

Also, it kind of seems like a good idea to put equations for the statistics into the little report.

##### 2026 Jul 29 (Wed)

Why is CvM not using expected order statistics? \* <https://claude.ai/share/13d39d37-a9ab-4331-9a33-253bf9d72e7b>
