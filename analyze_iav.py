# Implementación:

# 1. is_significant()


def is_significant(log2_fold_change, padj, lfc_threshold, padj_threshold):
    """Evalua si un gen cumple con los criterios de significancia."""

    if padj < padj_threshold and abs(log2_fold_change) >= lfc_threshold:
        return True

    return False


# ejemplo para probar is_significant():

print(is_significant(4.2, 0.0001, 1, 0.05))
print(is_significant(0.3, 0.0001, 1, 0.05))
print(is_significant(3.0, 0.8, 1, 0.05))


# 2. classify_gene()


def classify_gene(log2_fold_change):
    """Clasifica un gen significativo."""

    if log2_fold_change > 0:
        return "upregulated"

    return "downregulated"

# ejemplo para probar classify_gene():
print(classify_gene(4.2))
print(classify_gene(-3.0))


# 3. load_deseq2_results()

def load_deseq2_results(input_file):

    genes = []

    with open(input_file, "r") as file:

        next(file)

        for line in file:

            line = line.strip()

            if not line:
                continue

            columns = line.split("\t")

            if len(columns) < 7:
                continue

            gene = columns[0]

            try:
                log2_fold_change = float(columns[2])
                padj = float(columns[6])

            except ValueError:
                continue

            genes.append((gene, log2_fold_change, padj))

    return genes


# ejemplo para probar load_deseq2_results():
# print(load_deseq2_results("data/iav_deseq2_results.tsv"))


# 4. write_results()


def filter_genes(genes, lfc_threshold, padj_threshold):
    """Filtra genes significativos y los clasifica."""

    filtered_genes = []

    for gene, log2_fold_change, padj in genes:

        if is_significant(log2_fold_change, padj, lfc_threshold, padj_threshold):

            status = classify_gene(log2_fold_change)

            filtered_genes.append((gene, log2_fold_change, padj, status))

    return filtered_genes

# ejemplo para probar filter_genes():

genes = load_deseq2_results("data/iav_deseq2_results.tsv")

resultados = filter_genes(genes, 1, 0.05)

print(resultados[:10])


# 5. write_results()

def write_results(filtered_genes, output_file):
    """Escribe los genes significativos en un archivo TSV."""

    with open(output_file, "w") as file:

        file.write("gene\tlog2FoldChange\tpadj\tstatus\n")

        for gene, log2_fold_change, padj, status in filtered_genes:

            file.write(f"{gene}\t{log2_fold_change}\t{padj}\t{status}\n")

# ejemplo para probar write_results():

genes = load_deseq2_results("data/iav_deseq2_results.tsv")

filtered_genes = filter_genes(genes, 1, 0.05)

write_results(filtered_genes, "results/iav_significant_genes.tsv")

print("Archivo generado correctamente")
