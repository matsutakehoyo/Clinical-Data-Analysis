# Clinical-Data-Analysis

The repository contains details and Stan programs for the analysis of clinical data, accompanying the paper "".

## Analysis of FST

FST measurements consist subjective responses at various light intensities for different light colors. The stimulus may or may not be present, so that responses may be positive, negative, false negative, or false positive. We converted these responses to a dichotomous outcome (success or fail) depending on whether participants responded correctly. The outcome of the i-th trial $y_i$  (success = 1 or fail = 0) can be modeled with a Bernoulli distribution with probability to answer correctly  $theta$. 

$$y_i \sim Bernoulli(\theta_i)$$

We used logistic regression to estimate the effect of parameters on the probability of response $θ_i$. In addition to the conventional logistic regression, we implemented a ‘guessing’ parameter α_eye for robustness, as responses were sometimes seemingly random (Reference: Doing Bayesian Data Analysis (3rd Edition), Chapter 21- John K. Kruschke). The probability of responding correctly is therefor given by a completely random process (probability = 0.5), and the logistic regression part with $α_eye$ representing the fraction of ‘guessing’ (0≤α_eye≤1).

$$θ_i=\frac{1}/{2} \alpha_{pat}+(1-\alpha_{pat} )logistic(\mu_i )$$

For the ‘guessing’ parameter α we used a broad prior which gives values over 0.5 very low but non-zero probability.

$$α_eye \sim Beta(1,9)$$

For the logistic regression, we used light intensity ($x_i$), patient ($\beta_pat$), eye ($\beta_eye$), and light color ($\beta_clr$) as predictors. Since we were interested in characteristic responses to particular combinations of eye and light color, we included an interactions term between eye and light color ($\beta_eyeXclr$). This allowed us to estimate the overall trends for patient, eye, color, taking into account the hierarchical structure of the data. A sum-to-zero constraint ($\Sigma\beta_k=0$,for $k=pat,eye,clr,eyeXclr$) was imposed on these predictors, and posteriors shown as offsets from the overall mean ($\beta_0).

$$\mu_i=\beta_0+\beta_{pat,i}+\beta_{eye,i}+\beta_{clr,i}+\beta_{eyeXclr,i}+\beta x_i$$

For $\beta_pat$, $\beta_eye$, $\beta_eyeXclr$ we used full Bayesian inference with a $Gamma(1.64,0.32)$ for hyperpriors which has a mean of 2 and standard deviation of 4, covering all likely values in log odds scale.

β_k^ ~Normal(0,σ_k )  
σ_k~Gamma(1.64,0.32),     for k=pat,eye,eyeXclr 
As there are only four colors, we used an informative half-t prior for the effect of color 
β_clr~Student(3,0,1)

Analysis of chromatic pupillometry
Chromatic pupillometry measurements consist of five consecutive time series measurements of pupil diameter changes to light exposures aimed at stimulating rods, cones (two conditions), or melanopsin. For the melanopsin stimulation, pupil changes of the first repeat was clearly different from the rest of the successive measurements, and we therefore separated these measurements to mela1 (repeat 1) and mela2 (repeats 2-5). The time series measurements were condensed into three key regions: before (mean value before light exposure, 0 – 200 ms), peak (peak value after light exposure, 200 – 2000 ms), and after (mean value of steady state after light exposure, 3000-6000 ms). The peak region was only defined for rod/cone stimulation as we do not expect, nor observed, a transient response for melanopsin stimulation. Measurements sometimes contained regions where the standard deviation was zero. These regions with constant values were excluded from the analysis as they represent areas where pupil detection failed. Before, peak, and after values were analyzed using a multivariate regression, using multivariate student-t distribution for robustness. 
[x_(beore,i),x_(peak,i),x_(after,i) ]=t_ν (μ_(before,i),μ_(peak,i),μ_(after,i),Σ_clr)
Patient, eye, and color were considered as predictors, as well as the interaction between eye and color, with broad hyperpriors. A sum-to-zero constraint (Σβ=0,for pat,eye,clr,eyeXclr) was imposed on predictor, and posteriors are shown as offsets from the overall mean (β_0).
(μ_i ) ⃗=(β_0 ) ⃗+β ⃗_( pat,i)+β ⃗_(eye,i)+β ⃗_(clr,i)+β ⃗_( eyeXclr,i)   
β ⃗_k  ~ Normal(0,σ ⃗_k )
 σ_k  ~ Normal(0,10), for k=pat,eye,clr,exeXclr
As “peak” is not defined for mela1 and mela2, correlation coefficients were estimated by marginalizing the observed components. The correlation coefficients Ω were estimated using Cholesky factorization and LKJ(2) prior.
Σ_clr=diag(σ ⃗ )×Ω_clr× diag(σ ⃗ )         (σ ⃗=[σ_before,σ_peak,σ_after])
Ω_clr=L_clr L_clr'
L_clr~LKJ(2)
For the degrees of freedom ν we used the Gamma(2,0.1) prior recommended by the Stan development team.

Analysis of correlation 
Association between features (sex, age, EPT, logMAR, Retina thickness, FST and chromatic pupillometry measurements) was examined with a Bayesian counterpart of Pearson’s correlation test by estimating the correlation coefficient of a multivariate distribution. We implemented our model with a multivariate t-distribution, instead of a multi normal distribution for robustness against outliers. A latent state was assumed to take into account the uncertainty for measurements where mean and standard deviation were available (retina thickness, FST, and chromatic pupillometry). For FST and chromatic pupillometry data, the posterior estimates (mean and standard deviation) of the respective analyses were used as data. For missing data, i.e. missing components of the multivariate outcome, we modeled the marginal distribution of the component that is observed. Finally, for EPT measurements data contained censored values (right-censored at 2.5mA). These measurements were incorporated as parameters constrained to values in the censored range (>2.5mA). 
x ⃗_i  ~ t_ν (μ ⃗,Σ )
Σ=diag(σ ⃗ )×Ω× diag(σ ⃗ )
Ω=LL'
L~LKJ(2)
We performed a second correlation analysis similar to the above, but conditioning the correlation coefficients on visual acuity was performe to further explore the association between FST and retinal thickness. We separated patients by visual acuity into four categories: FZ (logMAR < 2.9), HB (logMar < 3.1), SL+(logMAR < 3.4), and SL- (logMAR ≥ 3.4), and estimated correlation coefficients within each group. 
x ⃗_i  ~ t_ν (μ ⃗,Σ_VA  )
Σ_VA=diag(σ ⃗ )×Ω_VA× diag(σ ⃗ )
Ω_VA=L_VA L_VA'
L~LKJ(2)
