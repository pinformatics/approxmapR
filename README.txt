NAME: README.TXT
DESCRIPTION: This file provides documentation for the R-package approxmapR.
Last Updated: 12/20/2021

===========================================
= Folders Descriptions (listed on GitHub) =
===========================================
data: Contains test datasets.
inst: Legend for reports.
man: contains the *.rd files for the package; these files are the how-to-use documentation for functions in the package.
R: The files which contain the R code that makes the package.
src:
tests: A few test cases which uses -testthat-.
vignettes: Contains examples on how to use approxmapR and helps one learn the package.


=======================
= Package DESCRIPTION =
=======================

::NOTE:: If using on a Windows computer, this package requires RTools to be installed on the machine since it uses C++ code for the computationally expensive portion. One can install RTools from https://cran.r-project.org/bin/windows/Rtools/rtools40.html.

Approxmap is an algorithm used for exploratory data analysis of sequential data which uses an approximate alignment approach developed by Dr. Hye-Chun Kum. This approach is used on longitudinal data when one is wanting to identify the underlying patterns.
ApproxmapR aims to provide a consistent and tidy API for using the algorithm in R. There are two vignettes which are encouraged to be reviewed before using the package; they can be found in the vignettes folder and are titled:

(1) approxmapR-getting-started.Rmd, and
(2) long_example.Rmd

Example data is provided while getting familiar with the package.


==============
= References =
==============

Kum, Hye-Chung, Pei, J., Wang, W., and Duncan, D. (2002). ApproxMAP: Approximate Mining of Consensus Sequential Patterns. Mining Sequential Patterns from Large Data Sets. Springer Series (28): The Kluwer International Series on Advances in Database Systems, pp. 138-160.
