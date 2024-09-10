import pandas as pd
from goatools import obo_parser
from goatools.associations import read_ncbi_gene2go
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy

# Load the GO DAG (ontology structure)
obodag = GODag("go-basic.obo")

# Load gene2go associations for human (taxid 9606)
geneid2gos_human = read_ncbi_gene2go("gene2go", taxids=[9606])

# Read genes of interest and background genes
with open("C://Users//tik105//Desktop//SOMSCAN//Figures//Figure 5 What do the chemicals and shRNAs do//Gene target list.txt") as f:
    genes_of_interest = [line.strip() for line in f]

with open("C://Users//tik105//Desktop//SOMSCAN//Figures//Figure 5 What do the chemicals and shRNAs do//Go background genes.txt") as f:
    background_genes = [line.strip() for line in f]

# Initialize GOEnrichmentStudy
go_enrichment_study = GOEnrichmentStudy(
    background_genes,  # List of background genes
    geneid2gos_human,  # geneid/GO associations
    obodag,            # Ontologies
    methods=['fdr_bh'] # Multiple test correction methods
)

# Run the enrichment analysis
results = go_enrichment_study.run_study(genes_of_interest)

# Print results
go_enrichment_study.print_results(results, min_overlap=1)


# Save results to a file
with open("go_enrichment_results_human.txt", "w") as fout:
    go_enrichment_study.print_results(results, min_overlap=1, pval=0.05, fout=fout)
