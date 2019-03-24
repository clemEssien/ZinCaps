import os

from glob import glob


DATA_DIR = "../../data/"
positive_pids = []
negative_pids = []
for files in glob(DATA_DIR+"*"):
    filename = os.path.split(files)[1]
    print(type(filename))
    file = open(DATA_DIR+filename,"r")
    for f in file:
        if "positive" in filename :
            positive_pids.append(f[0:6]);
        else:
            negative_pids.append(f[0:6])
    print("------------------------------")

print(len(negative_pids))
print(len(positive_pids))