# Pending allele decisions

Allele decisions recorded (DT, 2026-10-28). Rule: alt = the allele that
gnomAD v4 NFE and 1000 Genomes CEU and GBR all report at consistent frequency.
LongevityMap is supporting evidence only (mixed strands). All 12 blocks below
are `alt_status=resolved` in `data/processed/variants_47_input_v2.csv`; alleles
were not changed.

Sources:
- `analysis/panel/allele_audit_47.csv` — GRCh38 ref, Ensembl alt list
- `data/processed/allele_evidence_47.csv` — gnomAD NFE site alleles, LongevityMap allele
- `data/processed/variants_47_input_v2.csv` — `alt_status` and existing notes
- `data/processed/manual_resolution_alleles.csv` — header only (no rows)
- Ensembl REST release-116 variation `populations` (CEU and GBR AFs), matching the panel cache

`manual_resolution_alleles.csv` has no data rows.

---

## rs28366003 (MT2A) — alt_status=resolved

- GRCh38 ref base: A
- Ensembl alleles: C|G
- gnomAD NFE AF (site alleles): C=9.552458629170082e-06|G=0.05635305174759087
- CEU AF: A=0.94949494949495|G=0.0505050505050505
- GBR AF: G=0.0274725274725275|A=0.972527472527473
- LongevityMap allele: C
- v2 note: C=longevitymap|G=eqtl+gnomad_nfe
- Decision: alt G (DT, 2026-10-28). LongevityMap allele is on the opposite strand.

## rs1800774 (CETP) — alt_status=resolved

- GRCh38 ref base: C
- Ensembl alleles: A|G|T
- gnomAD NFE AF (site alleles): CTT=8.479505882233231e-07|T=0.358176615993499
- CEU AF: C=0.651515151515151|T=0.348484848484849
- GBR AF: C=0.620879120879121|T=0.379120879120879
- LongevityMap allele: T
- v2 note: T=longevitymap+gnomad_nfe
- Decision: alt T (DT, 2026-10-28).

## rs5882 (CETP) — alt_status=resolved

- GRCh38 ref base: G
- Ensembl alleles: A|C
- gnomAD NFE AF (site alleles): A=0.6847960989901594
- CEU AF: G=0.358585858585859|A=0.641414141414141
- GBR AF: G=0.296703296703297|A=0.703296703296703
- LongevityMap allele: (empty in allele_evidence_47.csv)
- v2 note: gnomAD v4 alternate alleles: A; Ensembl alternate alleles: A|C
- Decision: alt A (DT, 2026-10-28). Ref G is the minor allele in Europeans; the Ensembl C allele is absent from gnomAD NFE.

## rs41383 (NLRC5) — alt_status=resolved

- GRCh38 ref base: T
- Ensembl alleles: A|C
- gnomAD NFE AF (site alleles): C=0.5034861143798541
- CEU AF: T=0.545454545454545|C=0.454545454545455
- GBR AF: C=0.472527472527473|T=0.527472527472527
- LongevityMap allele: C
- v2 note: C=longevitymap+gnomad_nfe
- Decision: alt C (DT, 2026-10-28).

## rs9916344 (PITPNM3) — alt_status=resolved

- GRCh38 ref base: C
- Ensembl alleles: A|G|T
- gnomAD NFE AF (site alleles): G=0|T=0.5093323127649553
- CEU AF: C=0.47979797979798|T=0.52020202020202
- GBR AF: T=0.521978021978022|C=0.478021978021978
- LongevityMap allele: T
- v2 note: T=eqtl+longevitymap+gnomad_nfe
- Decision: alt T (DT, 2026-10-28).

## rs7207422 (ASIC2) — alt_status=resolved

- GRCh38 ref base: A
- Ensembl alleles: C|G
- gnomAD NFE AF (site alleles): G=0.2489481272251155
- CEU AF: G=0.227272727272727|A=0.772727272727273
- GBR AF: A=0.791208791208791|G=0.208791208791209
- LongevityMap allele: G
- v2 note: G=longevitymap+gnomad_nfe
- Decision: alt G (DT, 2026-10-28).

## rs1801318 (NDUFS1) — alt_status=resolved

- GRCh38 ref base: T
- Ensembl alleles: A|C
- gnomAD NFE AF (site alleles): C=0.308998458380533
- CEU AF: C=0.282828282828283|T=0.717171717171717
- GBR AF: T=0.642857142857143|C=0.357142857142857
- LongevityMap allele: (empty in allele_evidence_47.csv)
- v2 note: C=eqtl+gnomad_nfe
- Decision: alt C (DT, 2026-10-28).

## rs1981429 (SDC4) — alt_status=resolved

- GRCh38 ref base: G
- Ensembl alleles: C|T
- gnomAD NFE AF (site alleles): T=0.5274121838471277
- CEU AF: G=0.525252525252525|T=0.474747474747475
- GBR AF: G=0.406593406593407|T=0.593406593406593
- LongevityMap allele: T
- v2 note: T=eqtl+longevitymap+gnomad_nfe
- Decision: alt T (DT, 2026-10-28).

## rs139170 (PARVG) — alt_status=resolved

- GRCh38 ref base: C
- Ensembl alleles: G|T
- gnomAD NFE AF (site alleles): T=0.457037799676423
- CEU AF: C=0.580808080808081|T=0.419191919191919
- GBR AF: C=0.516483516483517|T=0.483516483516484
- LongevityMap allele: T
- v2 note: T=eqtl+longevitymap+gnomad_nfe
- Decision: alt T (DT, 2026-10-28).

## rs155979 (PCSK1) — alt_status=resolved

- GRCh38 ref base: G
- Ensembl alleles: A|C|T
- gnomAD NFE AF (site alleles): A=0|C=0.8582767439809507|T=0
- CEU AF: G=0.156565656565657|C=0.843434343434343
- GBR AF: G=0.142857142857143|C=0.857142857142857
- LongevityMap allele: C
- v2 note: C=longevitymap+gnomad_nfe
- Decision: alt C (DT, 2026-10-28). G/C palindrome resolved by frequency (MAF 0.14).

## rs524533 (TMEM151B) — alt_status=resolved

- GRCh38 ref base: C
- Ensembl alleles: A|G|T
- gnomAD NFE AF (site alleles): T=0.6948957965383257
- CEU AF: T=0.707070707070707|C=0.292929292929293
- GBR AF: T=0.653846153846154|C=0.346153846153846
- LongevityMap allele: T
- v2 note: T=longevitymap+gnomad_nfe
- Decision: alt T (DT, 2026-10-28).

## rs882696 (TUSC3) — alt_status=resolved

- GRCh38 ref base: T
- Ensembl alleles: A|C|G
- gnomAD NFE AF (site alleles): C=0.8738496932515337
- CEU AF: C=0.888888888888889|T=0.111111111111111
- GBR AF: T=0.142857142857143|C=0.857142857142857
- LongevityMap allele: A
- v2 note: A=longevitymap|C=eqtl+gnomad_nfe
- Decision: alt C (DT, 2026-10-28). LongevityMap allele is on the opposite strand.
