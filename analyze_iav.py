# Implementación:

# 1. is_significant()


def is_significant(log2_fold_change, padj, lfc_threshold, padj_threshold):
    """Evalua si un gen cumple con los criterios de significancia."""

    if padj < padj_threshold and abs(log2_fold_change) >= lfc_threshold:
        return True

    return False


# ejemplo para probar la funcion:

print(is_significant(4.2, 0.0001, 1, 0.05))
print(is_significant(0.3, 0.0001, 1, 0.05))
print(is_significant(3.0, 0.8, 1, 0.05))
