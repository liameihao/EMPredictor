from PromoterSVM.SVM.model import *
import sys
# Input number of features for each run
steps = []

# Name of the data
# Must be in the same format as the test.csv file:
# Label, Feature 1, Feature 2, ...
df_name = sys.argv[1]

# First, run the fselect algorithm to find the fscores for each feature
# fselect
fselect_main(df_name)
# Generate final files for fscores and features, is a dependency for post_fs_svm function
gen_finfscore(df_name)

# SVM after fselect, input the number of the features, the function itself will subset the original dataset into
# new dataset and run libsvm algorithm. The result is generated into the .result file.
for s in steps:
    print "number of features: ", s
    # SVM
    post_fs_svm(df_name, df_name+".finfscore", s)
