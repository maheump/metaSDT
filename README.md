metaSDT
=======

DESCRIPTION

metaSDT is a R function analyzing, with signal detection theory, your type II data. It works with confidence ratings (any type of scale with 1 step between two grades) and control levels (such as information seeking level or amount of exploration).

INPUTS

It only requires:
    - an "input_data" data.frame with answers labels in the first column and confidence ratings in the second one,
    - and confidence borns.
It can also take into account design (detection or choice), ponderation type, and homogeneity analysis, which are specified by default.

Answers have to be classified (answers labels) following classical SDT categories:
    1: Hit,
    2: False alarm,
    3: Miss,
    4: Correct rejection.
    
UTILIZATION

1) Open R
2) Set the working directory to the right place
3) Load the function
4) Load the example data
5) Launch the function

REFERENCES

Here are some usefull references to undestand what is going on here :
    Fleming, S. M., Weil, R. S., Nagy, Z., Dolan, R. J., & Rees, G. (2010). Relating introspective accuracy to individual differences in brain structure. Science, 329(5998), 1541–1543.
    Galvin, S. J., Podd, J. V., Drga, V., & Whitmore, J. (2003). Type 2 tasks in the theory of signal detectability: Discrimination between correct and incorrect decisions. Psychonomic Bulletin & Review, 10(4), 843–876.
    Kornbrot, D. E. (2006). Signal detection theory, the approach of choice: Model-based and distribution-free measures and evaluation. Perception & Psychophysics, 68(3), 393–414.
    Macmillan, N. A., & Creelman, C. D. (2004). Detection theory: A user's guide. Psychology press.
    Maniscalco, B., & Lau, H. (2012). A signal detection theoretic approach for estimating metacognitive sensitivity from confidence ratings. Consciousness and Cognition, 21(1), 422–430.
    Szczepanowski, R., & Pessoa, L. (2007). Fear perception: Can objective and subjective awareness measures be dissociated? Journal of Vision, 7(4), 1–17.