## Features
- Replace rate in each cluster
- Top replace genes in each cluster

## Output
- `{sample}.rep_in_tsne.txt` Replace rate in each cluster.
- `{sample}.rep_in_tsne_top10` Top 10 replace genes in each cluster.
## Arguments
`--tsne` tsne file from analysis step.

`--mat` matrix replacement file, from replacement step.

`--rep` cell replacement file, from replacement step.

`--mincell` turn-over in at least cells, default 5.

`--topgene` show top N genes,default 10.

`--outdir` Output diretory.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, celescope may output addtional file for debugging.

