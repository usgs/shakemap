.. _sec-verification-4:

****************************
Verification
****************************

The results produced by the ShakeMap **model** module are the product of
an interpolation scheme based on the statistics of multivariate
normal distributions (MVN). See :ref:`Worden et al. (2018) <worden2018>`
for a discussion of this approach. The mathematical complexity of 
the MVN process makes 
it difficult to ever fully verify the software against all possible 
inputs, or to even assert with certainty that any particular result is
objectively correct (at least once the inputs exceed some minimum 
level of complexity). Here, we discuss a set of simplified verification
tests that provide some support for the belief that the software is
producing correct results that are consistent with our hand calculations.
These tests
are not designed to fully test all of the features of the software --
that job is left to our unit tests and integration tests. Here we make
numerous simplifications in order to more easily interpret the results.

While the tests discussed in this section are one-dimensional (i.e.,
the results are computed for a line along which the sources are located), 
the computational process is agnostic to dimensionality and is only 
concerned
with the distances between locations. Again, our other testing considers
more complex models, and the results of those tests appear consistent
with the results presented here.

Various simplifying assumptions were made when producing these tests 
in order to better illuminate the behavior of the MVN process itself. 
In particular, the ground-motion prediction equation (GMPE) used
in these tests always returns a mean of 0 (in log space) for all locations, 
and reports a between-event standard deviation of 0.6 and a 
within-event standard deviation of 0.8 (making the total 
standard deviation a convenient 1.0). In addition, the 
cross-correlation function employed in these tests returns the product 
of the ratio of the
spectral periods (that is, ``Ts/Tl`` where ``Ts`` is the smaller period 
and ``Tl`` is the larger) and ``exp(-h/10)``, in which ``h`` is the 
separation distance. This model, while not the result of an empirical 
study, provides a smoother, more predictable behavior than other models
found in the literature and implemented in ShakeMap.

The verification tests may be run from the ShakeMap *bin* directory with 
the command **run_verification**. The command will run the tests and then
attempt to open a window displaying the plots. This last step might 
not work on all systems. The plots can be found in
*tests/data/eventdata/verification_test_XXXX/current/products* (where
"*XXXX*" is the number of the test).

Test 0001
====================

:num:`Figure #verification-test-one` shows the results of Test 0001. This
test places two observation points along a line. 
As discussed above, the GMPE evaluates to 0 (in log units) everywhere.  
Both observations in this test also have an amplitude of 0.0 (in log units), 
so the computed bias of the event is 0.
Thus, the conditional mean amplitude evaluates to 0 everywhere. The standard 
deviation is 0 at the location of the observations, and at great distances
from the observations it asymptotes to a value somewhat less than 1 (but
still greater than the GMPE's within-event standard deviaiton of 0.8).
This is because with only two observations, the considerable uncertainty
of the bias is applied to the within-event uncertainty.
These are the expected results, are consistent with our hand calculations,
and provide some confidence that the
MVN implementation is not introducing a bias or other anomalies.


.. _verification-test-one:

.. figure:: _static/verification_test_0001_PGA.*
   :width: 700
   :align: left

   Verification Test 0001. Two observations along a line have 
   amplitudes of 0.
   The black line shows the conditional mean, the blue lines
   show the conditional mean +/-- the conditional standard
   deviation (upper plot), and the red line shows the conditional
   standard deviation (lower plot).


Test 0002
====================

Test 0002 is shown in :num:`Figure #verification-test-two`. In this test,
one observation has an amplitude of +1.0, the other is --1.0. Because of
the offsetting observations, the bias is again 0. The figure shows that
the conditional amplitudes reach the expected value (+/-- 1.0) at the 
observation points, and approach 0 at distances far from the 
observations. As with Test 0001, the standard deviation is 0 at 
the observations and reaches a maximum somewhere between 0.8 and 1.0
at great distance from the observations.


.. _verification-test-two:

.. figure:: _static/verification_test_0002_PGA.*
   :width: 700
   :align: left

   Verification Test 0002. Two observations along a line have 
   amplitudes of +1.0 and --1.0.
   The black line shows the conditional mean, the blue lines
   show the conditional mean +/-- the conditional standard
   deviation (upper plot), and the red line shows the
   conditional standard deviation (lower plot).

Test 0003
====================

Verification Test 0003 has a single observation with an amplitude of +1.0
along a line (see :num:`Figure #verification-test-three`). 
If we apply equation 11 of :ref:`Worden et al. (2018) <worden2018>` we 
can determine the expected bias:

.. math::

    \mu_{\delta B_{i,m}} = 
    \frac{\bm{Z}^T_i \mathbf{\Sigma^{-1}}_{\mathbf{Y_2Y_2},i} 
    \bm{\zeta}_i }
    {\tau_{i, m}^{-2} + \bm{Z}^T_i \mathbf{\Sigma^{-1}}_{\mathbf{Y_2 Y_2},i}
    \bm{Z}_i},

In the bivariate case, this reduces to:

.. math::

    \mu_{\delta B_{i,m}} = 
    \frac{z \frac{1}{\phi^2} \zeta}  
    {\tau_{i, m}^{-2} + z \frac{1}{\phi^2} z}

In our case the GMPE mean is 0 and the observation is 1, giving 
:math:`\zeta = 1.0`. Because the observed IMT is the same as the
desired output IMT, the factors :math:`z` are 1.0. The within-event
standard deviation (:math:`\phi`) is 0.8, 
and the between-event standard deviation (:math:`\tau`) is 0.6. Thus
we have

.. math::

    \mu_{\delta B_{i,m}} = 
    \frac{1 \times \frac{1}{0.8^2} \times 1}  
    {\frac{1}{0.6^2} + 1 \times \frac{1}{0.8^2} \times 1}

or,


.. math::

    \mu_{\delta B_{i,m}} = 0.36

Thus, the bias is 0.36, as seen in :num:`Figure #verification-test-three` 
(solid black line) at distance from the observation.

The standard deviation of the bias is given by equation 12 of 
:ref:`Worden et al. (2018) <worden2018>`. The bivariate form is:

.. math::

    \sigma_{\delta B_{i,m}} = 
    \frac{1}  
    {\tau_{i, m}^{-2} + 1 \times \frac{1}{\phi^2} \times 1},

which, in our case, is

.. math::

    \sigma^2_{\delta B_{i,m}} = 
    \frac{1}  
    {\frac{1}{0.6^2} + \frac{1}{0.8^2}},

or,

.. math::

    \sigma_{\delta B_{i,m}} = 0.48.

Thus, the conditional standard deviation at great distance from 
an observation will be (as given by equation 13 of 
:ref:`Worden et al. (2018) <worden2018>`):

.. math::

    \sigma_{\delta W_{i,m}} = \sqrt{\phi^2_{i,m} + \sigma^2_{\delta B_{i,m}}},

or,

.. math::

    \sigma_{\delta W_{i,m}} = \sqrt{0.8^2 + 0.48^2} = 0.933,


as we see in :num:`Figure #verification-test-three` (red dashed line)
and can verify from the data produced by ShakeMap. 


.. _verification-test-three:

.. figure:: _static/verification_test_0003_PGA.*
   :width: 700
   :align: left

   Verification Test 0003. A single observation along a line with 
   an amplitude of +1.0.
   The black line shows the conditional mean, the blue lines
   show the conditional mean +/-- the conditional standard
   deviation (upper plot), and the red line shows the conditional
   standard deviation (lower plot).

Test 0004
====================

Test 0004 uses an identical set up to Test 0003, except there
are two observations (of amplitude +1.0) at the same location.
Because the observations are co-located and of the same period,
their correlation is 1.0. This means that they will have the
effect of a single observation. The result, illustrated in
:num:`Figure #verification-test-four` confirms this. Note that
:num:`Figure #verification-test-four` (which has two observations)
is identical to :num:`Figure #verification-test-three` (which
has only one observation).


.. _verification-test-four:

.. figure:: _static/verification_test_0004_PGA.*
   :width: 700
   :align: left

   Verification Test 0004. Two observations at the same 
   location along a line, both with 
   amplitudes of +1.0.
   The black line shows the conditional mean, the blue lines
   show the conditional mean +/-- the conditional standard
   deviation (upper plot), and the red line shows the conditional
   standard deviation (lower plot). Compare with 
   :num:`Figure #verification-test-three`.


Test 0004b
====================

Test 0004b uses an identical set up to Test 0004, except that the
two observations (of amplitude +1.0) have been separated by 1
degree of longitude. Thus, they are no longer highly correlated
and, consequently, the event bias has increased. 
The result is visualized in
:num:`Figure #verification-test-four-b` which may be compared with
:num:`Figure #verification-test-four`. Note that in Test 0004, the
conditional mean far from the observations was less than 0.5, 
whereas, in Test 4b, it is greater than 0.5; this consequence is
a result of the greater bias of Test 0004b. Similarly, the 
uncertainty at distance from the observations is slightly less
in Test 0004b than in Test 0004 because the two essentially
independent observations have reduced the uncertainty of the
bias, which lowers the overall uncertainty.


.. _verification-test-four-b:

.. figure:: _static/verification_test_0004b_PGA.*
   :width: 700
   :align: left

   Verification Test 0004b. Two observations at different
   locations along a line, both with amplitudes of +1.0.
   The black line shows the conditional mean, the blue lines
   show the conditional mean +/-- the conditional standard
   deviation (upper plot), and the red line shows the conditional
   standard deviation (lower plot). Compare with 
   :num:`Figure #verification-test-four`.
   

Test 0005
====================

Test 0005 also has two co-located observations (see Verification
Test 0004, above), but here they have
opposite amplitudes of +1.0 and --1.0. The result, shown in
:num:`figure #verification-test-five`, is that the conditional mean
and standard deviation behave as if there were only a single
observation with the mean amplitude of the two observations (i.e.,
0).


.. _verification-test-five:

.. figure:: _static/verification_test_0005_PGA.*
   :width: 700
   :align: left

   Verification Test 0005. Two observations at the same 
   location along a line, with amplitudes of +1.0 and --1.0.
   The black line shows the conditional mean, the blue lines
   show the conditional mean +/-- the conditional standard
   deviation (upper plot), and the red line shows the 
   conditional standard deviation (lower plot).


Test 0006
====================

:num:`Figure #verification-test-six` illustrates Verification Test 0006. 
Forty evenly-spaced observations, all with amplitudes of +1.0 are used. 
Most of the observations are to the left of center of the plot (and
extend some ways off the left edge of the plot). Here we note that 
the bias has moved significantly toward the mean of the data (as 
compared with a single observation as in 
:num:`Figure #verification-test-three`), and the conditional
standard deviation at distance has decreased toward the within-event
value of 0.8.


.. _verification-test-six:

.. figure:: _static/verification_test_0006_PGA.*
   :width: 700
   :align: left

   Verification Test 0006. Forty evenly-space observations along 
   a line, with amplitudes of +1.0 (note that the observations
   extend some distance off the left edge of the figure).
   The black line shows the conditional mean, the blue lines
   show the conditional mean +/-- the conditional standard
   deviation (upper plot), and the red line shows the conditional
   standard deviation (lower plot).

Test 0007
====================

Verification Test 0007 uses a single observation with an amplitude
of +1.0. The observation is of spectral acceleration (SA) at a 
period of 1.0 s. The conditional mean SA was 
computed for the location of the observation at a variety of 
periods ranging from 0.1 s to 10.0
s. A separate bias is computed for each period, and the
correlation between the observation and the amplitude being
computed decreases as the ratio of the two periods decreases,
thus the amplitude tends toward zero as the ratio of the periods
decreases. At periods far from the observation period, the 
bias approaches 0 and its standard deviation approaches the
between-event standard deviation, thus the conditional standard
deviation approaches the combined between-event and within-event
standard deviation (which, in our tests is 1.0).


.. _verification-test-seven:

.. figure:: _static/verification_test_0007_spectra_plot.*
   :width: 700
   :align: left

   Verification Test 0007. A single observation of spectral
   acceleration (with an amplitude of 1.0) at a period of
   1.0 s, and the conditional spectral acceleration
   at periods from 0.1 s to 10.0 s.
   The black line shows the conditional mean, the blue lines
   show the conditional mean +/-- the conditional standard
   deviation (upper plot), and the red line shows the 
   conditional standard deviation (lower plot).

Test 0008
====================

Verification Test 0008 demonstrates the effect of uncertainty in the
value at the observation point. 
:num:`Figure #verification-test-eight`
illustrates five separate cases to show
the effect of five values of additional standard deviation: 0.0, 0.75,
1.5, 3.0, and 6.0 on an observation with an amplitude of 1.0 (as in 
Test 0003).  As we did with Test 0003, we can compute the bias and 
the adjusted within-event standard deviation for each of the five cases.

.. math::

    \mu_{\delta B_{i,m}} = 
    \frac{z \frac{1}{\phi^2} \zeta}  
    {\tau_{i, m}^{-2} + z \frac{1}{\phi^2} z}

Here, however, the factors :math:`z` and the residual must be 
mutiplited by the "omega factors" to compensate for the 
uncertainty in the observation (see Equation 32 of 
:ref:`Worden et al. (2018) <worden2018>`):

.. math::

    \omega = \sqrt{\frac{\sigma^2}{\sigma^2 + \sigma^2_\epsilon}}

In our case, for the five values of :math:`\sigma_\epsilon`, and using
a within-event standard deviation of :math:`\sigma = 0.8`, 
:math:`\omega` is approximately 1.0, 0.73, 0.47, 0.26, and 0.13. As
we saw in 
Test 0003, the bias for the case of :math:`\sigma_\epsilon = 0` is
0.36. For the second case, where :math:`\sigma_\epsilon = 0.75` we
have:

.. math::

    \mu_{\delta B_{i,m}} = 
    \frac{0.73 \times \frac{1}{0.8^2} \times 0.73}  
    {\frac{1}{0.6^2} + 0.73 \times \frac{1}{0.8^2} \times 0.73}

or,

.. math::

    \mu_{\delta B_{i,m}} = 0.23

By application of the same process, we find the bias values for the
other cases. Thus, the bias values for the five cases are about 0.36,
0.23, 0.11, 0.037, and 0.01.
These values can be observed in the top plot in 
:num:`Figure #verification-test-eight`
at distance from the observation point.

Similarly the standard deviation of the bias,

.. math::

    \sigma_{\delta B_{i,m}} = 
    \frac{1}  
    {\tau_{i, m}^{-2} + z \frac{1}{\phi^2} z},

for each of the cases is approximately 0.48, 0.53, 0.57, 0.59,
and 0.60. Thus, using the formula:

.. math::

    \sigma_{\delta W_{i,m}} = \sqrt{\phi^2_{i,m} + \sigma^2_{\delta B_{i,m}}},

we get values for the revised within-event standard deviation of 0.93, 0.96,
0.98, 0.99, and 1.0. These values are illustrated in the bottom half of
:num:`Figure #verification-test-eight`
at distance from the observation point.

Equations 42 and 43 of :ref:`Worden et al. (2018) <worden2018>`
apply to the bivariate case where the observation is at the site we
wish to estimate. The mean and variance are given by:

.. math::

    \mu | z = \mu + \frac{\sigma^2}{\sigma^2 + \sigma^2_\epsilon}\left( z - \mu \right),

and

.. math::

    \sigma^2 | z = \frac{\sigma^2 \sigma^2_\epsilon}{\sigma^2 + \sigma^2_\epsilon}.

Note that in this case, we wish to use the revised within-event standard 
deviation for :math:`\sigma`, and the bias value for :math:`\mu` (since the mean
value returned by the GMPE is 0). Again, for the case of
:math:`\sigma_\epsilon = 0.75`, we have

.. math::

    \mu | z = 0.23 + \frac{0.96^2}{0.96^2 + 0.75^2}\left( 1.0 - 0.23 \right),

or

.. math::

    \mu | z = 0.71.

Then

.. math::

    \sigma | z = \sqrt{\frac{0.96^2 \times 0.75^2}{0.96^2 + 0.75^2}},

or

.. math::

    \sigma | z = 0.59

Thus, at 
the observation point we get mean values for each of our five cases of 1.0, 0.71,
0.38, 0.13, and 0.04. These values can be seen for each of the five lines in 
:num:`Figure #verification-test-eight`
at the location of the observation point. Similarly, the standard deviation 
at the observation
point for each case is 0.0, 0.59, 0.82, 0.94, and 0.99. These values also
correspond to what we see in the lower plot of 
:num:`Figure #verification-test-eight`.


.. _verification-test-eight:

.. figure:: _static/verification_test_0008_PGA.*
   :width: 700
   :align: left

   Verification Test 0008. Five separates runs of ShakeMap
   each using a single observation with an
   amplitude of +1.0, but with increasing uncertainty.
   The upper plot (solid lines) shows the conditional means,
   and the lower plot (dashed lines) shows the conditional standard
   deviations. The black lines should be identical to their 
   counterparts in :num:`Figure #verification-test-three`.

Test 0009
====================

Test 0009 (see 
:num:`Figure #verification-test-nine`) has five observations:
the central observation has an amplitude of 0.75, while the 
other four observations have amplitudes of 1.0. All five have 
a standard 
deviation of 0.2. The spacing of the higher amplitudes was 
chosen to exert a strong influence on the amplitude at the 
location of the central observation,
but despite that influence its conditional mean should approach 
the observational amplitude (0.75) from below, but not reach or 
exceed it.


.. _verification-test-nine:

.. figure:: _static/verification_test_0009_PGA.*
   :width: 700
   :align: left

   Verification Test 0009. Five observations: the central
   observation has an amplitude of 0.75, while the other
   four have amplitudes of 1.0. All five observations have
   a standard deviation of 0.2.
   The black line shows the conditional mean, the blue lines
   show the conditional mean +/-- the conditional standard
   deviation (upper plot), and the red line shows the 
   conditional standard deviation (lower plot).

Test 0010
====================

Like Test 0009, Test 0010 (see 
:num:`Figure #verification-test-ten`) has five observations:
the central observation has an amplitude of 0.75, while the 
other four observations have amplitudes of 1.0. All five have 
a standard 
deviation of 0.2. Here, though, the spacing of the higher 
amplitudes was chosen so that the conditional amplitude at 
the location of
the central observation would approach the assigned amplitude
from above. The amplitude should not (quite) reach the 
observational value (0.75), or go below it.


.. _verification-test-ten:

.. figure:: _static/verification_test_0010_PGA.*
   :width: 700
   :align: left

   Verification Test 0010. Five observations: the central
   observation has an amplitude of 0.75, while the other
   four have amplitudes of 1.0. All five observations have
   a standard deviation of 0.2.
   The black line shows the conditional mean, the blue lines
   show the conditional mean +/-- the conditional standard
   deviation (upper plot), and the red line shows the conditional
   standard deviation (lower plot). Compare with
   :num:`Figure #verification-test-nine`.

Test 0011
====================

The purpose of this test is to verify the orientation of the Vs30
grid and the generic amplification factors. The origin and magnitude
are those of the January 17, 1994, Northridge, California earthquake.
:num:`Figure #verification-test-eleven` is an image of 3.0 s PSA. It
shows that the coastline and 
other geographic features of the Vs30 map are in the proper orientation.
This test also uses two generic amplification files that cover the 
same geographic area: one file has values of 1.0 
for the northern half of the grid, and 0.0 for the southern half, while
the second file has values of 1.0 for the western half and 0.0 for 
the eastern half. Thus, the northwest quadrant has a combined 
amplification of 2.0, the northeast and southwest quadrants have
amplification factors of 1.0, and the southeast quadrant has an
amplification of 0.0. The figure demonstrates that the combined
amplifications are working correctly and are in the proper 
orientation.


.. _verification-test-eleven:

.. figure:: _static/verification_test_0011_PSA3p0.*
   :width: 700
   :align: left

   Verification Test 0011. 3 s PSA  map using the epicenter and magnitude
   of the January 17, 1994, Northridge, California earthquake. The
   coastline and other background features are the result of
   site amplification from the Vs30 file. The major north-south and
   east-west divisions are the result of generic amplification 
   factors.
