from mvpa2.tutorial_suite import *

#  every row vector in the data matrix becomes an observation or a sample in the dataset, and every column vector represents an individual variable or a feature. The concepts 
# of samples and features are essential for a dataset, hence we take a further, closer look.
data = [[  1,  1, -1],
         [  2,  0,  0],
         [  3,  1,  1],
         [  4,  0, -1]]
ds = Dataset(data)
print ds.shape

print "length of ds: %d" % len(ds)

print "Number of features: %d" % ds.nfeatures


