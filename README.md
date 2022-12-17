# calex

Deriviation and verification of adjusted calcium coefficients

This R shiny app has been designed as a tool to assist in the process
of determining suitable coefficients to use to estimated the albumin
adjusted calcium concentration. It assumes a linear relationship between
serum albumin and serum calcium concentration of the form calcium = 
slope x albumin + intercept. Using the slope so derived the adjusted
calcium concentration is equal to measured calcium concentration + 
(slope x (albumin set point - measured albumin concentration)).

