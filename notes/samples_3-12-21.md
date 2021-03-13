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

Yes that seems right. Looks like silencing was done via siRNA treatment

```
Treatment with siRNA was carried out 24 h after seeding 200,000 wildtype ...
```

Treated cells were used for looking at O fragments.

```
Cells with inactivated DNA Ligase 1 and 3 were used to map the distribution of Okazaki fragments.
```

Deleted everything in raw data and replaces with just the human mapping from
the above geo link.

Downloaded geo stuff to local PC and then uploaded to Crick server

```
sudo scp -i  ~/.ssh/lsanz_id_rsa /home/ethan/Downloads/GSE134224_RAW.tar \
lsanz@crick.cse.ucdavis.edu:/home/lsanz/eth/projects/GLOE-seq/rawdata
```



