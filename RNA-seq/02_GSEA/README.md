# Gene Set Enrichment Analysis (GSEA) and visualisation

## Running GSEA

The [GSEA](https://www.gsea-msigdb.org/gsea/index.jsp) tool was run with the following relevant parameters:



```
param	num	100
param	scoring\_scheme	weighted
param	norm	meandiv
param	mode	Max\_probe
param	include\_only\_symbols	true
param	set\_max	500
param	nperm	1000
param	order	descending
param	rnd\_seed	timestamp
param	set\_min	15
param	create\_svgs	false
param	sort	real
param	create\_gcts	false
param	help	false
param	save\_rnd\_lists	false
param	median	false
param	metric	Signal2Noise
param	make\_sets	true
param	rnd\_type	no\_balance
param	gui	false
param	permute	gene\_set
param	collapse	No\_Collapse
```

For all runs, the stu\_PGSC-ITAG-merged.gmt file was used (see ../input/) and experiment-specific cls files (see ../input/Monocultures/ and ../input/Interactions/).


See [GSEA User Guide](https://docs.gsea-msigdb.org/#GSEA/GSEA_User_Guide/) for instructions on running the tool. 


## Visualisation

Significant results were plotted using the `gseaFromStats` and `gseaPlot` functions from the [biokit](https://github.com/martingarridorc/biokit) package.
