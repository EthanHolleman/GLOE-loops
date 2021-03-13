# Figuring out sample stuff

Currently trying to figure out what samples are which and what all these
bigwig files represent. 

Downloaded all files from GLOE-seq paper GEO id [at this link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134225)
and got directory full of bigwig files which I am guessing are some
kind of representation of the break points but not really sure
how to interpret these yet.

Can basically just ignore the yeast stuff for now.

# Mapping of ligase compentent cells

I think this is what we really want

[GEO link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134224)

Here they did in human cells and with ligase and
ligase knockout human cells (?).

Yes that seems right.

```
Cells with inactivated DNA Ligase 1 and 3 were used to map the distribution of Okazaki fragments.
```

Deleted everything in raw data and replaces with just the human mapping from
the above geo link.

