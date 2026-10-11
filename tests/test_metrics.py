from bioaudit.metrics import classification_metrics


def test_metrics_handle_missing_class_in_predictions():
    result = classification_metrics(
        ["a", "a", "b", "c"],
        ["a", "a", "b", "b"],
        labels=["a", "b", "c"],
    )
    assert result["per_class"]["c"]["recall"] == 0.0
    assert "c" in result["missing_classes_in_pred"]
    assert result["macro_f1"] >= 0.0
    assert result["confusion_matrix"]["labels"] == ["a", "b", "c"]


def test_single_class_partition():
    result = classification_metrics(["a", "a"], ["a", "a"], labels=["a", "b"])
    assert result["macro_f1"] >= 0.0
    assert result["auroc"] is None
    assert result["missing_classes_in_true"] == ["b"]
