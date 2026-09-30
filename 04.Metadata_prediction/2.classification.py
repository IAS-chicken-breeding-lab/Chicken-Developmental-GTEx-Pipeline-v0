import os
import sys
import numpy as np
import pandas as pd
import classification_function as cf


wkdir = sys.argv[1]
group_method = sys.argv[2]
num_splits = int(sys.argv[3])
num_repeats = int(sys.argv[4])
feature_method = sys.argv[5]
num_features = int(sys.argv[6])
class_model = sys.argv[7]
class_label = sys.argv[8]
folder_path = sys.argv[9]
input_X_file = sys.argv[10]
input_y_file = sys.argv[11]
input_X_test_file = sys.argv[12]

input_path = os.path.join(wkdir, "data")
output_path = os.path.join(wkdir, folder_path)
os.makedirs(output_path, exist_ok=True)


def read_expression(path):
    df = pd.read_csv(path, sep=r"\s+", header=0, index_col=0)
    sample_ids = df.index.astype(str).tolist()
    genes = df.columns.astype(str).tolist()
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0)
    return np.log2(df.values + 1), sample_ids, genes


def read_labels(path, sample_ids):
    lab = pd.read_csv(path, sep="\t", header=0)
    if "label" not in lab.columns:
        raise ValueError("Label file must contain a 'label' column")

    if lab.shape[1] > 1:
        id_col = [c for c in lab.columns if c != "label"][0]
        lab[id_col] = lab[id_col].astype(str)
        y = lab.set_index(id_col).reindex(sample_ids)["label"]
        if y.isna().any():
            raise ValueError("Some expression samples are missing labels")
        return y.values

    if len(lab) != len(sample_ids):
        raise ValueError("Number of labels does not match number of samples")
    return lab["label"].values


X, sample_ids, gene_ids = read_expression(os.path.join(input_path, input_X_file))
y = read_labels(os.path.join(input_path, input_y_file), sample_ids)

prefix = f"{feature_method}.F{num_features}.{class_model}.L_{class_label}"
metrics_file = os.path.join(output_path, prefix + ".csv")
features_file = os.path.join(output_path, prefix + ".feature.csv")
predict_file = os.path.join(output_path, prefix + ".predict.csv")


if class_label != "Unknown":
    metric_rows, feature_rows, pred_rows = [], [], []
    splits = cf.split_data(X, y, method=group_method, num_splits=num_splits, num_repeats=num_repeats)

    for i, (train_idx, test_idx) in enumerate(splits):
        repeat = i // num_splits + 1
        fold = i % num_splits + 1

        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        X_train_sel, X_test_sel, idx, fs_time = cf.select_features(
            X_train, y_train, X_test, feature_method, num_features
        )
        y_pred, _, pred_time = cf.fit_predict(X_train_sel, y_train, X_test_sel, class_model)
        metrics = cf.evaluate(y_test, y_pred)

        metric_rows.append({"Repeat": repeat, "Fold": fold, **metrics,
                            "feature_select_time": fs_time, "predict_time": pred_time})
        feature_rows.extend({"Repeat": repeat, "Fold": fold, "gene_id": gene_ids[j]} for j in idx)
        pred_rows.extend({"Repeat": repeat, "Fold": fold, "ind_id": sample_ids[j],
                          "True_label": y_test[k], "Pred_label": y_pred[k]}
                         for k, j in enumerate(test_idx))

    pd.DataFrame(metric_rows).to_csv(metrics_file, index=False)
    pd.DataFrame(feature_rows).to_csv(features_file, index=False)
    pd.DataFrame(pred_rows).to_csv(predict_file, index=False)

else:
    if input_X_test_file in {"", "NA", "None"}:
        raise ValueError("input_X_test_file is required when class_label=Unknown")

    X_test, test_ids, test_genes = read_expression(os.path.join(input_path, input_X_test_file))
    if test_genes != gene_ids:
        raise ValueError("Training and test expression matrices must contain the same genes in the same order")

    X_train_sel, X_test_sel, idx, fs_time = cf.select_features(
        X, y, X_test, feature_method, num_features
    )
    y_pred, _, pred_time = cf.fit_predict(X_train_sel, y, X_test_sel, class_model)

    pd.DataFrame({
        "gene_id": [gene_ids[j] for j in idx]
    }).to_csv(features_file, index=False)

    pd.DataFrame({
        "ind_id": test_ids,
        "Pred_label": y_pred
    }).to_csv(predict_file, index=False)

    pd.DataFrame([{
        "Training_samples": len(sample_ids),
        "Test_samples": len(test_ids),
        "Selected_features": len(idx),
        "feature_select_time": fs_time,
        "predict_time": pred_time
    }]).to_csv(metrics_file, index=False)

print("Done:", output_path)
