from sklearn.feature_selection import mutual_info_classif
from sklearn.metrics import mutual_info_score
import numpy as np


a = np.array([1, 1, 1, 0, 0, 1, 0, 0, 0, 1])
b = np.array([1, 1, 1, 0, 0, 1, 0, 0, 0, 1])

print(mutual_info_classif(a.reshape(-1,1), b, discrete_features = True)) # mutual information of 0.69, expressed in nats
print(mutual_info_score(a,b))