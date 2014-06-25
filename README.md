metaSDT
=======

DESCRIPTION

metaSDT is a R function analyzing, with signal detection theory, your type II data. It works with confidence ratings (any type of scale with 1 step between two grades) and control levels (such as information seeking level or amount of exploration).

REQUIRED INPUTS

It only requires:
    - An "input_data" data frame with, at least, type I correction (1 vs. 0) in the first column and metacognitive score (confidence ratings, information seeking levels, etc.) in the second one.
    - Confidence borns ("min_conf" and "max_conf").
    
OPTIONAL INPUTS

It can also take into account:
   - Type I signal detection analysis if (i) answers are labelled following classical SDT categories (1: hit, 2: false alarm, 3: miss, 4: correct rejection) in a third column included in the input data frame and (ii) the corresponding argument is activated ("typeI = 1").
   - Design correction: detection ("present" vs. "absent": "design = 1") or choice ("first" vs. "second": "design = 2" by default).
   - Ponderation type: classical Aroc computation (1/2 + 1/4*ka + 1/4*kb: "coefficient = 1") or balanced Aroc computation (with k coefficients based on the (un)balance between low and high ratings: "coefficient = 2").
   - Output type: simply display plots and results ("output = 0"), save the results table in a variable you defined ("output = 1"), save the bias curve coordinates in a variable you defined ("output = 2") or save the ROC curve coordinates in a variable you defined ("output = 3").

REMARKS   
   
This function automatically compute the Aroc in both half of the trials before analyzing them at once. The first two Aroc are compared against each other in order to get an error and so provided an homogenity analysis of the metacognitive accuracy during the experiment.
    
UTILIZATION

1) Open R,
2) Set the working directory to the right place,
3) Load the function,
4) Load the example data,
5) Launch the function.

REFERENCES

Here are some usefull references to undestand what is going on here :
    Fleming, S. M., Weil, R. S., Nagy, Z., Dolan, R. J., & Rees, G. (2010). Relating introspective accuracy to individual differences in brain structure. Science, 329(5998), 1541–1543.
    Fleming, S. M., & Lau, H. C. (2014). How to measure metacognition. Frontiers in Human Neuroscience, 8, 443.
    Galvin, S. J., Podd, J. V., Drga, V., & Whitmore, J. (2003). Type 2 tasks in the theory of signal detectability: Discrimination between correct and incorrect decisions. Psychonomic Bulletin & Review, 10(4), 843–876.
    Kornbrot, D. E. (2006). Signal detection theory, the approach of choice: Model-based and distribution-free measures and evaluation. Perception & Psychophysics, 68(3), 393–414.
    Macmillan, N. A., & Creelman, C. D. (2004). Detection theory: A user's guide. Psychology press.
    Maniscalco, B., & Lau, H. (2012). A signal detection theoretic approach for estimating metacognitive sensitivity from confidence ratings. Consciousness and Cognition, 21(1), 422–430.
    Szczepanowski, R., & Pessoa, L. (2007). Fear perception: Can objective and subjective awareness measures be dissociated? Journal of Vision, 7(4), 1–17.