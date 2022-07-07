# Clinical-Data-Analysis

The repository contains details and Stan programs for the analysis of clinical data, accompanying the paper **Clinical Diagnosys**.

## Analysis of FST

FST measurements consist subjective responses at various light intensities for different light colors. The stimulus may or may not be present, so that responses may be positive, negative, false negative, or false positive. We converted these responses to a dichotomous outcome (success or fail) depending on whether participants responded correctly. The outcome of the i-th trial $y_i$, with success = 1 or fail = 0, can be modeled with a Bernoulli distribution with probability to answer correctly $　\theta_i　$. 

$$y_i \sim Bernoulli(\theta_i)$$

We used logistic regression to estimate the effect of parameters on the probability of response $θ_i$. In addition to the conventional logistic regression, we implemented a ‘guessing’ parameter α_eye for robustness, as responses were sometimes seemingly random (Reference: Doing Bayesian Data Analysis (3rd Edition), Chapter 21- John K. Kruschke). The probability of responding correctly is therefor given by a completely random process (probability = 0.5), and the logistic regression part with $α_eye$ representing the fraction of ‘guessing’ (0≤α_eye≤1).

$$θ_i=\frac{1}{2} \alpha_{pat}+(1-\alpha_{pat} )logistic(\mu_i )$$

For the ‘guessing’ parameter α we used a broad prior which gives values over 0.5 very low but non-zero probability.

$$α_eye \sim Beta(1,9)$$

For the logistic regression, we used light intensity ($x_i$), patient ($\beta_{pat}$), eye ($\beta_{eye}$), and light color ($\beta_{clr}$) as predictors. Since we were interested in characteristic responses to particular combinations of eye and light color, we included an interactions term between eye and light color ($\beta_{eyeXclr}$). This allowed us to estimate the overall trends for patient, eye, color, taking into account the hierarchical structure of the data. A sum-to-zero constraint ($\Sigma\beta_k=0$,for $k=pat,eye,clr,eyeXclr$) was imposed on these predictors, and posteriors shown as offsets from the overall mean ($\beta_0).

$$\mu_i=\beta_0+\beta_{pat,i}+\beta_{eye,i}+\beta_{clr,i}+\beta_{eyeXclr,i}+\beta x_i$$

For $\beta_{pat}$, $\beta_{eye}$, $\beta_{eyeXclr}$ we used full Bayesian inference with a $Gamma(1.64,0.32)$ for hyperpriors which has a mean of 2 and standard deviation of 4, covering all likely values in log odds scale.

$$
\beta \sim Normal(0, \sigma_k) \\
\sigma_k \sim Gamma(1.64, 0.32), for k = pat, eye, eyeXclr \\
$$

As there are only four colors, we used an informative half-t prior for the effect of color 

$$\beta_{clr}~Student(3,0,1)

## Analysis of chromatic pupillometry

Chromatic pupillometry measurements consist of five consecutive time series measurements of pupil diameter changes to light exposures aimed at stimulating rods, cones (two conditions), or melanopsin. For the melanopsin stimulation, pupil changes of the first repeat was clearly different from the rest of the successive measurements, and we therefore separated these measurements to mela1 (repeat 1) and mela2 (repeats 2-5). The time series measurements were condensed into three key regions: before (mean value before light exposure, 0 – 200 ms), peak (peak value after light exposure, 200 – 2000 ms), and after (mean value of steady state after light exposure, 3000-6000 ms). The peak region was only defined for rod/cone stimulation as we do not expect, nor observed, a transient response for melanopsin stimulation. Measurements sometimes contained regions where the standard deviation was zero. These regions with constant values were excluded from the analysis as they represent areas where pupil detection failed. Before, peak, and after values were analyzed using a multivariate regression, using multivariate student-t distribution for robustness. 

$$[x_{beore,i},x_{peak,i},x_{after,i} ] = t_{\nu} (μ_{before,i},μ_{peak,i},μ_{after,i}, \Sigma_{clr})$$

Patient, eye, and color were considered as predictors, as well as the interaction between eye and color, with broad hyperpriors. A sum-to-zero constraint ($\Sigma\beta=0$,for pat,eye,clr,eyeXclr) was imposed on predictor, and posteriors are shown as offsets from the overall mean ($\beta_0$).

$$\vec{\mu_i} = \vec{\beta_0} + \vec{\beta_{pat,i}} + \vec{\beta_{eye,i}} + \vec{\beta_{clr,i}} + \vec{\beta_{eyeXclr,i}}$$
$$\beta_k \sim Normal(0, \sigma_k)$$  
$$\sigma_k \sim Normal(0,10), for k=pat,eye,clr,exeXclr $$

As “peak” is not defined for mela1 and mela2, correlation coefficients were estimated by marginalizing the observed components. The correlation coefficients $\Omega$ were estimated using Cholesky factorization and LKJ(2) prior.
$$\Sigma_{clr} = diag(\vec{\delta}) \times \Omega_{clr} \times diag(\vec{\delta}), (\vec{\sigma} = [\sigma_{before}, \sigma_{peak}, \sigma_{after}]))$$
$$\Omega_{clr} = L_{clr} L_{clr}'$$
$$L_{clr} \sim LKJ(2)$$

For the degrees of freedom ν we used the $Gamma(2,0.1)$ prior recommended by the Stan development team.

## Analysis of correlation 

Association between features (sex, age, EPT, logMAR, Retina thickness, FST and chromatic pupillometry measurements) was examined with a Bayesian counterpart of Pearson’s correlation test by estimating the correlation coefficient of a multivariate distribution. We implemented our model with a multivariate t-distribution, instead of a multi normal distribution for robustness against outliers. A latent state was assumed to take into account the uncertainty for measurements where mean and standard deviation were available (retina thickness, FST, and chromatic pupillometry). For FST and chromatic pupillometry data, the posterior estimates (mean and standard deviation) of the respective analyses were used as data. For missing data, i.e. missing components of the multivariate outcome, we modeled the marginal distribution of the component that is observed. Finally, for EPT measurements data contained censored values (right-censored at 2.5mA). These measurements were incorporated as parameters constrained to values in the censored range (>2.5mA). 

$$\vec{x_i} \sim t_{nu}(\vec{\mu}, \Sigma) $$
$$\Sigma = diag(\vec{\delta}) \times \Omega \times diag(\vec{\delta})$$
$$\Omega = LL'$$
$$L \sim LKJ(2)$$

We also performed a second correlation analysis similar to the above, but conditioning the correlation coefficients on visual acuity was performe to further explore the association between FST and retinal thickness. We separated patients by visual acuity into four categories: CF (logMAR < 2.9), HM (logMar < 3.1), LP(logMAR < 3.4), and NLP (logMAR ≥ 3.4), and estimated correlation coefficients within each group. 

$$\vec{x_i} \sim t_{nu}(\vec{\mu}, \Sigma_{VA}) $$
$$\Sigma_{VA} = diag(\vec{\delta}) \times \Omega_{VA} \times diag(\vec{\delta})
$$\Omega_{VA} = L_{VA}L_{VA}'$$

## Quadratic Regression

A regression analysis was performed to investigate the relationship between FST values (for example Blue vs Red). As Blue FST values generally exhibited the lowes threshold values, we used Blue as a reference, and compared how Green, White, and Red FST values deviate from the Blue values. While FST values exhibit an overall linear trend, the Blue and Red FST values could not be satisfactorily approximated  with a line, as Red values deviated more from Blue values at lower thresholds, and linear regression resluted in physiologically unrealistic values at high threhold regions (i.e. Red FST value much lower than Blue values at high thresholds ~40dB). We therefore used a quadratic regression model. This results in a curve in which Blue FST values never exceed Red FST values which is more physiologically consistent. The predictive posterior distribution of the theshold value (value at which probaility = 0.1) from the FST analysis were used as inputs and outputs for the regression. We attempted to implement a mesurement error model, which takes into acount both unncertainty in inputs $x$ (Blue FST) and outputs $y$ (White, Green, Red FST) values, however mcmc sampling did not converge when both measurement errors were presnt. We therefor implemented the measurement error only for $x$ (Blue FST values), assuming a latent state $x_{lat}$. 

$$  x \sim Normal(x_{lat}, \sigma_x) $$

Where $x$ and $\sigma_x$ are the mean and standard deviation of the posterior predictie distribution of the Blue FST value, and $x_lat$ is a latant variable. Thus the FST values for Green, White, and Red measurements ($y_i$) are given by the coefficients for the quadratic regression ($b_{0,clr}$, $b_{1,clr}$, and $b_{2,clr}$.

$$ \mu_i = b_{0,clr} + b_{1,clr} * x_{lat,i} + b_{2,clr} * (x_{lat,i})^2 $$
$$ y_i \sim Normal(\mu_i, \sigma_{clr}) $$


