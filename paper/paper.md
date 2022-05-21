---
title: 'chyper: Functions for Conditional Hypergeometric Distributions'
tags:
  - R
  - conditional hypergeometric distributions
  - classification
  - method evaluation
authors:
  - name: William Nickols
    orcid: 0000-0001-8214-9746
    affiliation: 1
affiliations:
 - name: William Nickols, Harvard College
   index: 1
date: 20 May 2022
bibliography: paper.bib

---

# Summary

Machine learning classification methods often attempt to assign labels to a new data point from a set of possible labels.  When multiple classification tools assign different sets of labels to the same data point and the correct classification is unknown, it can be unclear whether the overlapping labels in the assignments are a reflection of the data point's true properties or simply a convergence by chance.  To answer this question, `chyper` is an R[@rmanual] package for working with conditional hypergeometric distributions, distributions describing the number of overlapping labels due to chance when $n$ sets of labels are assigned by $n$ classification tools.  In this framework, classification tool $i$ has a population of labels in its database with size $n_i$ not shared among all tools, a population of labels with size $s$ shared in the databases of all tools, and a sample size of $m_i$ representing how many labels it assigned to a data point.  The `chyper` package is avilable from CRAN [here](https://cran.r-project.org/package=chyper), and it can be installed in R with the command `install.packages("chyper")`.

# Statement of need

The package `chyper` is an R package for implementing conditional hypergeometric distributions to aid in comparing classification tools.  While conditional or hierarchical hypergeometric models have been described in bioinformatics literature [@kim_how_2018], an implementation in R has not yet been produced.  Like other R packages built for distributions, `chyper` includes a probability mass function, a cumulative distribution function, a quantile function, and a random number generator.  In addition, it implements functions to give the mean of a conditional hypergeometric distribution, a p-value for observing a particular number of overlaps, a maximum likelihood estimator (MLE) for an unknown database overlap size ($s$), an MLE and method of moments (MOM) estimator for an unknown database non-overlap size ($n_i$), and an MLE and MOM estimator for an unknown sample size ($m_i$).

This package was designed to be used in comparing classification tools for problems where each data point can receive multiple labels and each classification tool assigns potentially many labels to each data point.  For example, in microbiome data, a single metagenomic sample can have many taxa assigned to it by tools with databases containing thousands of taxa, but the assignments will likely differ by tool, and the databases will also differ by tool [@sun_challenges_2021].  Alternatively, image classification often involves methods assigning multiple tags to an image, but different tools might assign different tags, and those different tools might differ in which tags are contained in their databases [@CZERNIAWSKI2020103131].  This package provides a mathematically sound way to quantify the probability that labels assigned to some data by multiple methods are due to some actual feature of the data rather than chance.

# Mathematical Details

Consider two overlapping populations with $n_1$ items unique to population 1, $n_2$ items unique to population 2, and $s$ itmems in both populations (for a total population 1 size of $n_1+s$ and total population 2 size of $n_2+s$).  Take a random sample of $m_1$ items from population 1 (all combinations equally likely) and a random sample of $m_2$ from population 2 (all combinations equally likely).  Each sample is internally without replacement (i.e. the sample from population 1 can't include the same item twice), but the same item can be chosen in both the sample from population 1 and the sample from population 2.  Let $X_1$ denote the number of items sampled from population 1 that are in the intersection.  Then, $X_1\sim\textrm{HGeom}(s,n_1,m_1)$.  Then, let $X_2$ denote the number of items sampled from the second population that were already sampled in the first population's sample.  Then, conditional on $X_1$, the number of overlaps $X_2$ is distributed as $X_2\sim\textrm{HGeom}(X_1,n_2+s-X_1,m_2)$.  By the law of total probability, $$P(X_2=k)=\sum_{x_1}P(X_2=k|X_1=x_1)P(X_1=x_1)$$$$=\sum_{x_1=k}^{\textrm{min}(m_1,s)}\frac{{{x_1}\choose{k}}{{n_2+s-x_1}\choose{m_2-k}}}{{{n_2+s}\choose{m_2}}}\cdot\frac{{{s}\choose{x_1}}{{n_1}\choose{m_1-x_1}}}{{{n_1+s}\choose{m_1}}}$$

## Conditional hypergeometric extension to arbitrarily many populations
This solution can be extended to the case of $k$ overlapping populations with $n_i$ items in population $i$ not in the intersection of all the populations and $s$ items in the intersection of all the populations.  Note that some items counted towards $n_i$ might also be counted towards $n_{j}$ for $i\neq j$; they just cannot be counted towards all other $n_j$ or they would be counted towards $s$.  As before, let $m_i$ items be sampled from each population without replacement internally but with replacement between the sampling of each population.  Let $X_k$ be the number of items chosen in all $k$ samples.  Then, let $X_1$ be the number of items from $s$ in $m_1$, $X_2$ the number of people from $X_1$ in $m_2$, $X_3$ be the number of people from $X_2$ in $m_3$, and so on.  Then, extending the law of total probability from earlier,

$$\begin{split}
P(X_2=x_2) =& \sum_{x_1}P(X_2=x_2|X_1=x_1)P(X_1=x_1) \\
P(X_3=x_3) =& \sum_{x_2}P(X_3=x_3|X_2=x_2)P(X_2=x_2)=\\
&\sum_{x_2}P(X_3=x_3|X_2=x_2)\Big[\sum_{x_1}P(X_2=x_2|X_1=x_1)P(X_1=x_1)\Big]\\
...\\
P(X_k=k) =& \sum_{x_{k-1}}P(X_k=k|X_{k-1}=x_{k-1})\Big[\sum_{x_{k-2}}P(X_{k-1}=x_{k-1}|X_{k-2}=x_{k-2})\cdot\Big[...\cdot\\
&\Big[\sum_{x_1}P(X_2=x_2|X_1=x_1)P(X_1=x_1)\Big]...\Big]\Big]\\
=& \sum_{x_{k-1}=k}^{\textrm{min}(m_{k-1},s)}\frac{{{x_{k-1}}\choose{k}}{{n_k+s-x_{k-1}}\choose{m_k-k}}}{{{n_{k}+s}\choose{m_k}}}\cdot\Big[\sum_{x_{k-2}=k}^{\textrm{min}(m_{k-2},s)}\frac{{{x_{k-2}}\choose{x_{k-1}}}{{n_{k-1}+s-x_{k-2}}\choose{m_{k-1}-x_{k-1}}}}{{{n_{k-1}+s}\choose{m_{k-1}}}}\cdot\Big[...\cdot\\
&\Big[\sum_{x_1=k}^{\textrm{min}(m_1,s)}\frac{{{x_1}\choose{k}}{{n_2+s-x_1}\choose{m_2-k}}}{{{n_2+s}\choose{m_2}}}\cdot\frac{{{s}\choose{x_1}}{{n_1}\choose{m_1-x_1}}}{{{n_1+s}\choose{m_1}}}\Big]...\Big]\Big]
\end{split}$$

## Computational optimizations
Two dynamic programming optimizations speed this computation significantly.  First, calculating this PMF can be optimized by storing all the $P(X_i = x_i)$ at each level because they are the same regardless of any level above them (i.e. $P(X_1 = x_1)$ does not depend on the value of $k$ in $P(X_2=k)$).  Therefore, by storing $P(X_1=x_1)$ for $x_1$ from $0$ to $\textrm{min}(m_1,s)$, $P(X_2=x_2)$ can be computed without needing to recompute anything for the $P(X_1=x_1)$ level.  This can be repeated with the $P(X_2=x_2)$ level versus the $P(X_3=x_3)$ level and so on.  

Second, the calculation can be optimized by storing the $P(X_{i} = x_{i}|X_{i-1} = x_{i-1})$ and updating from $P(X_i = x_i|X_{i-1} = x_{i-1})$ to $P(X_i = x_i+1|X_{i-1} = x_{i-1})$.  For example, $$P(X_2=x_2)=\sum_{x_1=x_2}^{\textrm{min}(m_1,s)}\frac{{{x_1}\choose{x_2}}{{n_2+s-x_1}\choose{m_2-x_2}}}{{{n_2+s}\choose{m_2}}}\cdot\frac{{{s}\choose{x_1}}{{n_1}\choose{m_1-x_1}}}{{{n_1+s}\choose{m_1}}}$$ can be updated to $$P(X_2=x_2+1)=\sum_{x_1=x_2+1}^{\textrm{min}(m_1,s)}\frac{{{x_1}\choose{x_2+1}}{{n_2+s-x_1}\choose{m_2-(x_2+1)}}}{{{n_2+s}\choose{m_2}}}\cdot\frac{{{s}\choose{x_1}}{{n_1}\choose{m_1-x_1}}}{{{n_1+s}\choose{m_1}}}$$ by storing the $\frac{{{s}\choose{x_1}}{{n_1}\choose{m_1-x_1}}}{{{n_1+s}\choose{m_1}}}$ term from $P(X_1=x_1)$ and shifting it as necessary to keep proper indexing.  

In the equation
$$c\cdot\frac{{{x_1}\choose{x_2}}{{n_2+s-x_1}\choose{m_2-x_2}}}{{{n_2+s}\choose{m_2}}} = \frac{{{x_1}\choose{x_2+1}}{{n_2+s-x_1}\choose{m_2-(x_2+1)}}}{{{n_2+s}\choose{m_2}}}$$  
$$c=\frac{(x_1-x_2)(m_2-x_2)}{(x_2+1)(n_2+s-x_1-(m_2-(x_2+1)))}$$  

This holds for any level; that is, multiplying $P(X_i=x_i|X_{i-1}=x_{i-1})$ by $$\frac{(x_{i-1}-x_i)(m_i-x_i)}{(x_i+1)(n_i+s-x_{i-1}-(m_i-(x_i+1)))}$$ gives $P(X_i=x_i+1|X_{i-1}=x_{i-1})$.  Storing $P(X_i=x_i|X_{i-1}=x_{i-1})$ and $P(X_{i-1}=x_{i-1})$ and updating them element-wise in the sum rather than recomputing the hypergeometric distribution dramatically speeds up the computation.  For more than two groups, this entire process is performed recursively (i.e. $P(X_1=x_1)$ is used to calculate $P(X_2=x_2)$, which is used to calculate $P(X_3=x_3)$ and so on).

# Acknowledgements

I’d like to thank Srihari Ganesh and Skyler Wu for their help in reviewing and improving the two-population case; Max Li for his help in optimizing the PMF calculation; and Meghan Short and Kelsey Thompson from the Huttenhower Lab at the Harvard T.H. Chan School of Public Health and the Broad Institute for their help in reviewing the math behind the arbitrary-population-number extension.

# References
