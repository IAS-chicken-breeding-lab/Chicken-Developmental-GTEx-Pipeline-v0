import time
from sklearn.model_selection import RepeatedKFold, RepeatedStratifiedKFold
from sklearn.feature_selection import SelectKBest, SelectPercentile, f_classif, chi2, mutual_info_classif, RFE, RFECV, SelectFromModel
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, AdaBoostClassifier, BaggingClassifier, ExtraTreesClassifier, VotingClassifier
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score


def split_data(X, y, method="stratified_kfold", num_splits=5, num_repeats=20, random_state=1000):
    if method == "k_fold":
        cv = RepeatedKFold(n_splits=num_splits, n_repeats=num_repeats, random_state=random_state)
        return cv.split(X)
    if method == "stratified_kfold":
        cv = RepeatedStratifiedKFold(n_splits=num_splits, n_repeats=num_repeats, random_state=random_state)
        return cv.split(X, y)
    raise ValueError("group_method must be 'k_fold' or 'stratified_kfold'")


def select_features(X_train, y_train, X_test, method="f_classif", num_features=10):
    start = time.time()

    if method == "f_classif":
        selector = SelectKBest(f_classif, k=num_features)
    elif method == "chi2":
        selector = SelectKBest(chi2, k=num_features)
    elif method == "mutual_info_classif":
        selector = SelectKBest(mutual_info_classif, k=num_features)
    elif method == "select_from_model":
        selector = SelectFromModel(LogisticRegression(max_iter=5000), max_features=num_features)
    elif method == "rfe":
        selector = RFE(LogisticRegression(max_iter=5000), n_features_to_select=num_features)
    elif method == "rfecv":
        selector = RFECV(LogisticRegression(max_iter=5000), step=1, cv=5)
    elif method == "percentile":
        selector = SelectPercentile(f_classif, percentile=num_features)
    else:
        raise ValueError("Unsupported feature selection method")

    X_train_sel = selector.fit_transform(X_train, y_train)
    X_test_sel = selector.transform(X_test)
    feature_idx = selector.get_support(indices=True)
    return X_train_sel, X_test_sel, feature_idx, time.time() - start


def build_model(model="logistic_regression", random_state=1000):
    if model == "logistic_regression":
        clf = LogisticRegression(max_iter=5000, n_jobs=-1)
    elif model == "svc":
        clf = SVC()
    elif model == "decision_tree":
        clf = DecisionTreeClassifier(random_state=random_state)
    elif model == "k_neighbors":
        clf = KNeighborsClassifier(n_jobs=-1)
    elif model == "gaussian_nb":
        clf = GaussianNB()
    elif model == "mlp":
        clf = MLPClassifier(max_iter=2000, random_state=random_state)
    elif model == "random_forest":
        clf = RandomForestClassifier(n_jobs=-1, random_state=random_state)
    elif model == "gradient_boosting":
        clf = GradientBoostingClassifier(random_state=random_state)
    elif model == "adaboost":
        clf = AdaBoostClassifier(random_state=random_state)
    elif model == "bagging":
        clf = BaggingClassifier(n_jobs=-1, random_state=random_state)
    elif model == "extra_trees":
        clf = ExtraTreesClassifier(n_jobs=-1, random_state=random_state)
    elif model == "voting":
        clf = VotingClassifier([
            ("lr", LogisticRegression(max_iter=5000)),
            ("rf", RandomForestClassifier(n_jobs=-1, random_state=random_state)),
            ("gb", GradientBoostingClassifier(random_state=random_state))
        ])
    else:
        raise ValueError("Unsupported classification model")

    return make_pipeline(StandardScaler(), clf)


def fit_predict(X_train, y_train, X_test, model="logistic_regression"):
    clf = build_model(model)
    start = time.time()
    clf.fit(X_train, y_train)
    pred = clf.predict(X_test)
    return pred, clf, time.time() - start


def evaluate(y_true, y_pred):
    return {
        "Accuracy": accuracy_score(y_true, y_pred),
        "Precision": precision_score(y_true, y_pred, average="weighted", zero_division=0),
        "Recall": recall_score(y_true, y_pred, average="weighted", zero_division=0),
        "F1": f1_score(y_true, y_pred, average="weighted", zero_division=0)
    }
