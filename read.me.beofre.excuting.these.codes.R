############
## Update ##
############

# 1), I found newer version of "rEDM" using different functions to excute the CCM and MDR s-map,
# 2), and some users seems no solutions to access the 1.2.3 version of "rEDM" package (i.e. the version I used in Nature communications 2023)
# below, I provide the codes to download the 1.2.3 version of "rEDM"
# if you got new issues, free to let me know.

packageurl <- "https://cran.r-project.org/src/contrib/Archive/rEDM/rEDM_1.2.3.tar.gz"
install.packages(packageurl, repos=NULL, type="source")

packageVersion("rEDM") #check the package version if it is 1.2.3
library("rEDM")# load the package, version 1.2.3