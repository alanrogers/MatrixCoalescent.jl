# MatrixCoalescent

This package calculates the probability that there are *k* ancestral
lineages at the ancient end of an epoch of population history, given
that there are *n* descendant lineages at the recent end. It also
calculates the expected length of the interval during which there are
*k* lineages, where 1 <= k <= n. The argument is based on a
continuous-time Markov chain described in appendix I of
[Tavare 1984](http://www.damtp.cam.ac.uk/user/st321/CV_&_Publications_files/STpapers-pdf/T84a.pdf),
[Griffiths and Tavare 1998](https://www.tandfonline.com/doi/abs/10.1080/15326349808807471),
and
[Wooding and Rogers 2002](https://www.genetics.org/content/161/4/1641.short).
The chain begins with a modern sample of *n* haploid lineages, sampled
at the recent end of the epoch. As we trace the ancestry of this
sample into the past, the original sample of *n* lineages falls to
*n-1*, then *n-2*, and so on until only a single lineage is left, or
we reach the end of the epoch.

This Markov chain is well known, but it is not often used because
accurate calculations are difficult with samples of even modest size.
However, it is possible to factor the calculations into two steps, one
of which can be done in exact arithmetic, and only needs to be done
once at the beginning of the computer program. Numerical error arises
only in the second step, and is modest up to samples of size 35. 
The current code cannot cope with samples larger than this, because of
integer overflow.

Within an epoch, the population has constant haploid size *2N*.  We
measure time backwards from the recent end of the epoch in units of
*2N* generations. On this scale, time is *v = t/2N*, where *t* is time
in generations. Calculations involve the type `MatCoal`, which can be
constructed using any type of floating-point number. Thus, one can do
high-precision work using `MatCoal{BigFloat}`, or use
`MatCoal{Float64}` to emphasize speed.  The package exports several
functions:

* `MatCoal` constructs an object of type `MatCoal`.

* `eigenvals!` calculates eigenvalues for use by `project!` and
`interval_lengths!`.

* `project!` calculates the probabilities of *k=2,3,...,n* lineages at
  the ancient end of the epoch, given that there are *n* lineages at
  the recent end. The probability that *k=1* is obtained by
  subtraction.

* `interval_lengths!` calculates the expected lengths of the
  subintervals during which there are *k=2,3,...,n* lineages. The
  expected length of the subinterval within which *k=1* is obtained by
  subtraction.

