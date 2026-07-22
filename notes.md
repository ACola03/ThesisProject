
##### 2026 Jul 08 (Wed)

Adam to send link to a nicer version of provisional notes

We have many uniformity statistics:
* Histogram variance (requires boxing)
* variance of spaces (analytically nice)
	* start at 0, end at 1 for n+1 spaces
* naive location
* weighted location

The location statistics seem terrible! This is a big surprise and we should probably try to understand it at least for a while. In the meantime, let's go forward with the spacing statistic.

##### 2026 Jul 03 (Fri)

Order-statistic estimators. Currently liking [the WLSE from here](https://arxiv.org/pdf/2107.09316): 

----------------------------------------------------------------------

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

- and then for LRT

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
