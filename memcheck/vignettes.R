# Code to be run with
#   R -d "valgrind --tool=memcheck --leak-check=full" --vanilla < tests/thisfile.R
# First build and install the package.
library('ape')
devtools::build_vignettes()
