# Shiny clonetree
Interactive visualization of clonal tree in R/Shiny

To start, use
```{bash}
docker-compose up
```

Expected input data:

RDS file with a Seuarat object
```
./shiny_seurat/apps/clonetree.app/data/seurat.rds
```

Lineage tree in newick format (parentheses and commas). Entities should match the cell names from the Seurat object
```
./shiny_seurat/apps/clonetree.app/data/tree.treefile
```
