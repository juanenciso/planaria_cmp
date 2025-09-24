# planaria_cmp — Marker families across planarian species

Pipeline to identify and compare germline/gene-silencing marker families
(**Vasa/DDX4, PL10/DDX3, Piwi, Argonaute/Ago, Nanos**) between **Dugesia** and
**Schmidtea**, and to generate summary tables and publication-ready figures.

## Repo layout

planaria_cmp/
├─ scripts/ # Bash/R helpers (optional)
├─ refs/
│ ├─ dug/transcripts.fa.gz
│ └─ smed/transcripts.fa.gz
├─ quants/
│ ├─ dug/<ACC>/quant.sf
│ └─ smed/<ACC>/quant.sf
├─ db/ # UniProt marker FASTA + DIAMOND DB
├─ hits/ # DIAMOND tabular results
├─ markers/ # Collapsed tables + per-family FASTA + comparisons
├─ figures/ # PNG/PDF figures
├─ environment.yml
└─ README.md

> **Note:** Only commit compact, derived files (TSV, PNG/PDF). Raw reads, large
FASTA/indices, and caches should be `.gitignore`-d.

## 1) Environment

```bash
micromamba env create -f environment.yml
micromamba activate planaria-cmp

This installs: diamond, seqkit, and R (with ggplot2, dplyr, readr, tidyr, scales, forcats, patchwork).
2) Input expectations
	•	Transcripts (per species):
	◦	refs/dug/transcripts.fa.gz
	◦	refs/smed/transcripts.fa.gz
	•	Salmon/Kallisto quantification:
	◦	quants/dug/<ACC>/quant.sf
	◦	quants/smed/<ACC>/quant.sf
Pick one representative accession per species (used to fetch TPM per transcript).
3) Build the UniProt marker database
mkdir -p db

# Collect reviewed UniProt entries for each family keyword
: > db/markers.fa
for q in \
  'gene_exact:DDX4' \
  'gene_exact:DDX3X' 'gene_exact:DDX3Y' 'gene_exact:PL10' \
  'gene_exact:PIWIL1' 'gene_exact:PIWIL2' 'gene_exact:PIWIL3' 'gene_exact:PIWIL4' \
  'gene_exact:AGO1' 'gene_exact:AGO2' 'gene_exact:AGO3' 'gene_exact:AGO4' \
  'gene_exact:NANOS1' 'gene_exact:NANOS2' 'gene_exact:NANOS3'
do
  echo "↓ $q"
  curl -sL "https://rest.uniprot.org/uniprotkb/search?query=${q}+AND+reviewed:true&format=fasta" >> db/markers.fa
done

# Quick sanity check (FASTA headers)
grep -c '^>' db/markers.fa

# Build DIAMOND DB
diamond makedb --in db/markers.fa -d db/markers.dmnd

4) DIAMOND searches

mkdir -p hits

# Dugesia
diamond blastx \
  -q refs/dug/transcripts.fa.gz \
  -d db/markers.dmnd \
  --max-target-seqs 50 -e 1e-5 -k 50 \
  --outfmt 6 qseqid sseqid pident length evalue bitscore stitle \
  > hits/dug_vs_markers.tsv

# Schmidtea
diamond blastx \
  -q refs/smed/transcripts.fa.gz \
  -d db/markers.dmnd \
  --max-target-seqs 50 -e 1e-5 -k 50 \
  --outfmt 6 qseqid sseqid pident length evalue bitscore stitle \
  > hits/smed_vs_markers.tsv

5) Collapse to best hit per family and add TPM
Replace <ACC_SM> and <ACC_DU> with your accessions.
accSM=<ACC_SM>   # e.g., SRR32973373
accDU=<ACC_DU>   # e.g., SRR27692727

# Helper: assign family from sseqid/title
fam_awklib='
function fam_of(sid, st, a,b){
  a=tolower(sid); b=tolower(st);
  if (a ~ /ddx4/  || b ~ /vasa|ddx4/)               return "Vasa/DDX4";
  if (a ~ /ddx3l|ddx3|pl10/ || b ~ /ddx3|pl10/)     return "PL10/DDX3";
  if (a ~ /piwil|piwl|piwi/ || b ~ /piwi/)          return "Piwi";
  if (a ~ /ago[0-9]?/ || b ~ /argonaute|ago[0-9]?/) return "Argonaute/Ago";
  if (a ~ /nanos|nano[0-9]?/ || b ~ /nanos|nano/)   return "Nanos";
  return "";
}'

mkdir -p markers

# Best-by-family (Schmidtea)
awk -vFS='\t' -vOFS='\t' "$fam_awklib
{
  q=$1; sid=$2; ev=$5+0; bs=$6+0; st=$7;
  fam=fam_of(sid, st); if (fam=='') next;
  key=q OFS fam;
  if (!(key in seen) || ev<bestE[key] || (ev==bestE[key] && bs>bestB[key])) {
    seen[key]=1; bestE[key]=ev; bestB[key]=bs; row[key]=$0; famName[key]=fam;
  }
}
END{ for (k in row) print row[k], famName[k]; }" \
  hits/smed_vs_markers.tsv > hits/smed_best_by_family.tsv

# Best-by-family (Dugesia)
awk -vFS='\t' -vOFS='\t' "$fam_awklib
{ q=$1; sid=$2; ev=$5+0; bs=$6+0; st=$7;
  fam=fam_of(sid, st); if (fam=='') next;
  key=q OFS fam;
  if (!(key in seen) || ev<bestE[key] || (ev==bestE[key] && bs>bestB[key])) {
    seen[key]=1; bestE[key]=ev; bestB[key]=bs; row[key]=$0; famName[key]=fam;
  }
}
END{ for (k in row) print row[k], famName[k]; }" \
  hits/dug_vs_markers.tsv > hits/dug_best_by_family.tsv

# Join TPM from quant.sf
awk -vFS='\t' -vOFS='\t' -v qsf="quants/smed/${accSM}/quant.sf" '
BEGIN{ while((getline<qsf)>0){ if(NR==1)continue; tpm[$1]=$5 }
       print "transcript_id","TPM","family","stitle" }
{ id=$1; st=$7; fam=$8; print id, ((id in tpm)?tpm[id]:0), fam, st }' \
  hits/smed_best_by_family.tsv | (read h; echo "$h"; sort -k2,2nr) \
  > markers/smed_markers_TPM.tsv

awk -vFS='\t' -vOFS='\t' -v qsf="quants/dug/${accDU}/quant.sf" '
BEGIN{ while((getline<qsf)>0){ if(NR==1)continue; tpm[$1]=$5 }
       print "transcript_id","TPM","family","stitle" }
{ id=$1; st=$7; fam=$8; print id, ((id in tpm)?tpm[id]:0), fam, st }' \
  hits/dug_best_by_family.tsv | (read h; echo "$h"; sort -k2,2nr) \
  > markers/dug_markers_TPM.tsv

6) Per-family FASTA (optional)

# IDs
awk -F'\t' 'NR>1{print > ("markers/" tolower($3) ".ids")}' markers/smed_markers_TPM.tsv
awk -F'\t' 'NR>1{print > ("markers/dug_" tolower($3) ".ids")}' markers/dug_markers_TPM.tsv

# FASTA
seqkit grep -f markers/vasa/ddx4.ids  refs/smed/transcripts.fa.gz > markers/smed_vasa.fa  || true
# (repeat for other families and Dugesia)

7) Summary tables for plotting

# Long table: TPM sum per family/species
{
  echo -e "species\tfamily\tTPM_sum\tn"
  awk -F'\t' 'NR>1{s+=$2; c[$3]+=$2; n[$3]++} END{for(k in c) print "Schmidtea",k,c[k],n[k]}' markers/smed_markers_TPM.tsv
  awk -F'\t' 'NR>1{c[$3]+=$2; n[$3]++} END{for(k in c) print "Dugesia",k,c[k],n[k]}'       markers/dug_markers_TPM.tsv
} | tr ' ' '\t' > markers/family_TPM_long.tsv

# Comparison table
awk -F'\t' '
NR==FNR && FNR>1{a[$2]=$3; an[$2]=$4; next}
FNR>1{b[$2]=$3; bn[$2]=$4}
END{
  print "family","smed_TPM","smed_n","dug_TPM","dug_n","log2fc"
  for (k in a){
    st=a[k]+0; dt=(k in b)?b[k]+0:0; sn=an[k]+0; dn=(k in bn)?bn[k]+0:0;
    l=log((st+1e-6)/(dt+1e-6))/log(2);
    print k,st,sn,dt,dn,l
  }
}' OFS='\t' markers/family_TPM_long.tsv markers/family_TPM_long.tsv \
| awk 'NR==1||$1=="Vasa/DDX4"||$1=="PL10/DDX3"||$1=="Piwi"||$1=="Argonaute/Ago"||$1=="Nanos"' \
> markers/compare_family_TPM.tsv

8) Plots
Open R and run:

# scripts/plot_markers.R equivalent (inline)
need <- c("readr","dplyr","tidyr","ggplot2","scales","forcats","patchwork")
inst <- need[!need %in% rownames(installed.packages())]
if (length(inst)>0) install.packages(inst, repos="https://cloud.r-project.org")
library(readr); library(dplyr); library(ggplot2); library(scales)
library(forcats); library(patchwork)

dir.create("figures", showWarnings = FALSE)
long <- readr::read_tsv("markers/family_TPM_long.tsv", show_col_types = FALSE)
comp <- readr::read_tsv("markers/compare_family_TPM.tsv", show_col_types = FALSE)

fam_levels <- c("Vasa/DDX4","PL10/DDX3","Piwi","Argonaute/Ago","Nanos")
pal <- c("Vasa/DDX4"="#4E79A7","PL10/DDX3"="#F28E2B","Piwi"="#59A14F",
         "Argonaute/Ago"="#E15759","Nanos"="#B07AA1")
lab_fun <- if (utils::packageVersion("scales") >= "1.2.0")
  scales::label_number(scale_cut = scales::cut_short_scale()) else scales::label_number_si()

long$family <- factor(long$family, levels=fam_levels)

p_abs_en <- ggplot(long, aes(family, TPM_sum, fill=family)) +
  geom_col(width=.75, color="grey20") +
  facet_wrap(~ species, nrow=1) +
  scale_fill_manual(values=pal, guide="none") +
  scale_y_continuous(labels=lab_fun) +
  labs(title="Total TPM by marker family", x=NULL, y="TPM (sum)") +
  theme_minimal(12) + theme(panel.grid.major.x = element_blank(),
                            axis.text.x = element_text(angle=20, hjust=1),
                            strip.text = element_text(face="bold"))

long_pct <- long %>% group_by(species) %>% mutate(TPM_pct = 100*TPM_sum/sum(TPM_sum)) %>% ungroup()
p_pct_en <- ggplot(long_pct, aes(family, TPM_pct, fill=family)) +
  geom_col(width=.75, color="grey20") +
  facet_wrap(~ species, nrow=1) +
  scale_fill_manual(values=pal, guide="none") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(title="Within-species composition by family", x=NULL, y="% of total TPM") +
  theme_minimal(12) + theme(panel.grid.major.x = element_blank(),
                            axis.text.x = element_text(angle=20, hjust=1),
                            strip.text = element_text(face="bold"))

comp2 <- comp %>% mutate(family=factor(family, levels=fam_levels))
p_fc_en <- ggplot(comp2, aes(y = forcats::fct_rev(family), x = log2fc, color=family)) +
  geom_vline(xintercept=0, linetype="dashed", color="grey50") +
  geom_segment(aes(x=0, xend=log2fc, yend=forcats::fct_rev(family)), linewidth=1) +
  geom_point(size=3, stroke=0) +
  scale_color_manual(values=pal, guide="none") +
  scale_x_continuous(breaks=scales::pretty_breaks()) +
  labs(title="log2FC (Schmidtea / Dugesia) by family",
       x="log2 fold-change (TPMsum Smed / TPMsum Dug)", y=NULL) +
  theme_minimal(12) + theme(panel.grid.major.y = element_blank())

# Save individual figures
ggsave("figures/family_TPM_sum.png", p_abs_en, width=9, height=4.2, dpi=300)
ggsave("figures/family_TPM_pct.png", p_pct_en, width=9, height=4.2, dpi=300)
ggsave("figures/family_log2FC_lollipop.png", p_fc_en, width=7.2, height=4.2, dpi=300)

# Composite
combo <- (p_abs_en / p_pct_en / p_fc_en) +
  patchwork::plot_annotation(
    title = "Germline/gene-silencing marker families across planarian species",
    tag_levels = "A",
    theme = theme(plot.title = element_text(face="bold", size=16))
  )
ggsave("figures/markers_family_composite_en.png", combo, width=10, height=12, dpi=300)
ggsave("figures/markers_family_composite_en.pdf", combo, width=10, height=12)

9) Troubleshooting

	•	diamond: command not found → activate the environment.
	•	Error opening file db/markers.dmnd → run diamond makedb in step 3.
	•	[WARN] 0 patterns loaded from file (seqkit) → empty .ids file; check family names in TSV.
	•	R: label_number_si() defunct → the script already switches to label_number(scale_cut=...).
License

---

### `environment.yml`

```yaml
name: planaria-cmp
channels:
  - conda-forge
  - bioconda
dependencies:
  # core CLI tools
  - diamond >=2.1
  - seqkit >=2.5
  - curl
  - pigz
  - parallel

  # R stack
  - r-base >=4.3
  - r-ggplot2
  - r-dplyr
  - r-tidyr
  - r-readr
  - r-scales
  - r-forcats
  - r-patchwork

  # optional: tidyverse meta (comment out if not needed)
  # - r-tidyverse






