# Clinical-Data-Analysis

When $a \ne 0$, there are two solutions to $(ax^2 + bx + c = 0)$ and they are 
$$ x = {-b \pm \sqrt{b^2-4ac} \over 2a} $$

The repository contains details and Stan programs for the analysis of clinical data, accompanying the paper **Clinical Diagnosys**.

## Analysis of FST

FST measurements consist subjective responses at various light intensities for different light colors. The stimulus may or may not be present, so that responses may be positive, negative, false negative, or false positive. We converted these responses to a dichotomous outcome (success or fail) depending on whether participants responded correctly. The outcome of the i-th trial $y_i$  (success = 1 or fail = 0) can be modeled with a Bernoulli distribution with probability to answer correctly. 



