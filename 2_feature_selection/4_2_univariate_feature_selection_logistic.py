import pickle
from feature_engine.selection import SelectBySingleFeaturePerformance
from sklearn.linear_model import LogisticRegression
import pandas as pd
from imblearn.over_sampling import SMOTE
from collections import Counter
import numpy as np

# Define the path to the .pkl file
input_path = "/home/sr933/rcc_data/scRCC/combined_data/RCC_data_dict.pkl"

# Load the dictionary from the .pkl file
with open(input_path, "rb") as f:
    data_dict = pickle.load(f)

# Access the contents of the dictionary
X_combined = data_dict["X"].T
y_labels = data_dict["y"]
gene_list = data_dict["Genes"]

# Set all elements >= 1 to 1
y_labels[y_labels >= 1] = 1


##Do SMOTE

def data_smote(x,y):
    oversample = SMOTE()
    x_new, y_new = oversample.fit_resample(x, y)

    counter_old = Counter(y)
    print(counter_old)
    counter = Counter(y_new)
    print(counter)

    #Make features binary
    x_new=np.rint(x_new)
    return x_new, y_new

X_train, y_train= data_smote(X_combined, y_labels)

print("Starting feature selection")

# Use simple logistic regression as the estimator
sel = SelectBySingleFeaturePerformance(
    variables=None,
    estimator=LogisticRegression(solver="liblinear", max_iter=1000),
    scoring="roc_auc",
    cv=5,
    threshold=0.5
)

# Fit to find predictive features
sel.fit(X_train, y_train)

# Extract feature performance
feat = sel.feature_performance_

# Convert to DataFrame and save as CSV
feat1 = pd.DataFrame.from_dict(feat, orient='index', columns=["AUC"])
feat1.to_csv("/home/sr933/rcc/data/Logistic_Features.csv")

print("Feature selection results saved to Logistic_Features.csv.")


