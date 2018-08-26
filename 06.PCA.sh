### This file contains command line for PCA using eigensoft (https://www.hsph.harvard.edu/alkes-price/software/)
module add eigensoft/6.0.1
/path/tp/eigensoft/bin/smartpca -p fly.par > fly.log
perl /path/tp/eigensoft/bin/ploteig2  -i fly.evec -c 1:2 -x -o fly.pca12.xtxt -p CS45:CS65:DS45:DS65:HS45:HS65:KS45:KS65:SS45:SS65:UC4:UC45:UC65
perl /path/tp/eigensoft/bin/ploteig2  -i fly.evec -c 3:4 -x -o fly.pca34.xtxt -p CS45:CS65:DS45:DS65:HS45:HS65:KS45:KS65:SS45:SS65:UC4:UC45:UC65

