#!/bin/csh -f
# Run a couple tests of the normal point processing. 
# Check the files in data/analysis against those in data/analysis/ref.
bin/ldb_crd s25y11d016t0231\#103
echo comparing normal point files:
diff data/analysis/s25y11d016t0231\#103.npt data/analysis/ref/s25y11d016t0231\#103.npt
echo comparing full rate files:
diff data/analysis/s25y11d016t0231\#103.frd data/analysis/ref/s25y11d016t0231\#103.frd
bin/ldb_crd s25y11d042t0234\#103
diff data/analysis/s25y11d042t0234\#103.npt data/analysis/ref/s25y11d042t0234\#103.npt
diff data/analysis/s25y11d042t0234\#103.frd data/analysis/ref/s25y11d042t0234\#103.frd

