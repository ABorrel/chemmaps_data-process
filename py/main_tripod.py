import DBrequest
import tripod


pr_tripod = "/mnt/c/Users/aborr/research/sandbox/Tox21_assays/tripod_11-2020/"


c_tripod = tripod.tripod(pr_tripod)

l_col = c_tripod.get_colnames()
print("==Col to add in DB==")
print(l_col)

c_tripod.pushInDB("tox21_tripod")