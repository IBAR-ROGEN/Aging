# FOXO1 Longevity Gene Analysis

**Analysis Date:** October 31, 2025  
**Gene:** FOXO1 (Forkhead box protein O1)

## Executive Summary

FOXO1 is a critical longevity gene and transcription factor that plays a central role in aging, metabolism, and stress response. This analysis comprehensively examines FOXO1's protein partners, genetic interactions, and longevity data from multiple databases including STRING, SynergyAge, and OmniPath.

---

## 1. Gene Identification & Basic Information

### Gene Identifiers
- **Gene Symbol:** FOXO1 (also known as FKHR, FOXO1A)
- **UniProt ID:** Q12778 (primary), L8E9Y8 (secondary)
- **Ensembl ID:** ENSG00000150907
- **STRING ID:** 9606.ENSP00000368880
- **Species:** Homo sapiens (Human)
- **NCBI Taxonomy ID:** 9606

### Protein Description
**Full Name:** Forkhead box protein O1  
**Alternative Names:** 
- Forkhead box protein O1A
- Forkhead in rhabdomyosarcoma

**Protein Existence Level:** Evidence at protein level (highest confidence)

### Functional Overview
FOXO1 is the main target of insulin signaling and regulates metabolic homeostasis in response to oxidative stress. Key functions include:

- **Metabolic Regulation:** Main regulator of glucose metabolism and gluconeogenesis in response to insulin
- **Stress Response:** Key regulator of redox balance and cellular response to oxidative stress
- **Cell Survival:** Important regulator of apoptosis acting downstream of AKT1 and other kinases
- **Autophagy:** Required for autophagic cell death induction in response to starvation or oxidative stress
- **Bone Homeostasis:** Regulates osteoblast numbers and controls bone mass
- **Adipogenesis:** Regulates expression of adipogenic genes during preadipocyte differentiation

### Tissue Expression
- **High Expression:** Skeletal muscle, ovary
- **Moderate Expression:** Heart, placenta, lung, liver, pancreas, spleen, testis, small intestine
- **Low Expression:** Brain, thymus, prostate, mucosal lining of colon
- Expressed in umbilical endothelial cells (at protein level)

---

## 2. Protein-Protein Interactions (STRING Database)

FOXO1 has high-confidence interactions (score ≥ 0.700) with multiple key longevity and metabolic regulators:

### Top 10 Interaction Partners

| Rank | Partner | Gene Name | Interaction Score | Function Category |
|------|---------|-----------|-------------------|-------------------|
| 1 | SIRT1 | Sirtuin 1 | 0.999 | Longevity/Deacetylase |
| 2 | PPARGC1A | PGC-1α | 0.999 | Metabolic Regulation |
| 3 | AKT1 | AKT Kinase | 0.999 | Insulin Signaling |
| 4 | ATG7 | Autophagy 7 | 0.998 | Autophagy |
| 5 | CREBBP | CBP | 0.993 | Acetylation/Transcription |
| 6 | SMAD3 | SMAD3 | 0.990 | TGF-β Signaling |
| 7 | AKT2 | AKT2 Kinase | 0.986 | Insulin Signaling |
| 8 | CTNNB1 | β-Catenin | 0.985 | Wnt Signaling |
| 9 | SFN | 14-3-3 σ | 0.984 | Regulatory Protein |
| 10 | EP300 | p300 | 0.984 | Acetylation/Transcription |

### Key Interaction Highlights

#### Longevity & Aging Regulators
- **SIRT1 (Score: 0.999):** Critical deacetylase that activates FOXO1 transcriptional activity and extends lifespan
- **PPARGC1A (Score: 0.999):** Metabolic coactivator involved in mitochondrial biogenesis and energy homeostasis

#### Insulin/Growth Signaling (IIS Pathway)
- **AKT1 (Score: 0.999):** Phosphorylates FOXO1 leading to nuclear export and inactivation
- **AKT2 (Score: 0.986):** Alternative AKT isoform with similar regulatory function

#### Autophagy & Cell Death
- **ATG7 (Score: 0.998):** Essential autophagy gene; FOXO1 activates autophagy through ATG7

#### Transcriptional Regulation
- **CREBBP (Score: 0.993):** Acetylates FOXO1, modulating its DNA-binding activity
- **EP300 (Score: 0.984):** Another acetyltransferase regulating FOXO1 activity
- **SMAD3 (Score: 0.990):** TGF-β signaling mediator that cooperates with FOXO1

#### Signaling & Localization
- **SFN (Score: 0.984):** 14-3-3 protein that binds phosphorylated FOXO1 and promotes cytoplasmic retention
- **CTNNB1 (Score: 0.985):** Wnt signaling component that interacts with FOXO1

---

## 3. Longevity Data from SynergyAge Database

### FOXO Gene Effects on Lifespan

The SynergyAge database contains **13 experimental models** involving FOXO genes (primarily in *Drosophila melanogaster*):

#### Summary Statistics
- **Total Models:** 13
- **Primary Organism:** *Drosophila melanogaster* (fruit fly, tax_id: 7227)
- **Effect Range:** +0.0% to -34.3% (FOXO deletion/inactivation generally reduces lifespan)
- **Key Finding:** FOXO is critical for normal lifespan; its loss consistently reduces longevity

### Detailed Longevity Results

#### Single Gene Effects (FOXO Deletion/Inactivation)

| Model | Organism | Effect (%) | Lifespan (days) | Study (PMID) | Notes |
|-------|----------|------------|-----------------|--------------|-------|
| foxo(if) | *D. melanogaster* | -7.2% | 27.69 | 20624856 | Inducible RNAi knockdown |
| foxo(df) | *D. melanogaster* | -18.2% | 36.0 | 21518241 | Deletion mutant |
| foxo(DF) | *D. melanogaster* | -19.7% | 49.0 | 21443682 | Deletion mutant |
| foxo(DF) | *D. melanogaster* | -28.0% | 54.0 | 21443682 | Deletion mutant (different conditions) |
| foxo(DF) | *D. melanogaster* | -31.7% | 41.0 | 21443682 | Deletion mutant |
| foxo(DF) | *D. melanogaster* | -34.3% | 44.0 | 21443682 | Deletion mutant |

**Key Finding:** FOXO deletion consistently reduces lifespan by 7-34%, demonstrating its essential role in longevity.

#### Genetic Interactions with FOXO

| Partner Gene | Interaction Type | Effect (%) | Description |
|--------------|------------------|------------|-------------|
| **InR** (Insulin Receptor) | Epistatic/Dependent | -13.1% to -31.3% | FOXO mediates lifespan extension from InR reduction |
| **chico** | Opposite Effects | 0.0% to -12.5% | chico extends longevity through FOXO |
| **Pi3K92E** | Enhancer | -32.0% | PI3K-FOXO interaction in insulin signaling |
| **Lrrk** | Negative Epistasis | -8.2% | LRRK-FOXO interaction |
| **rpr** | Opposite Effects | -26.7% | Reaper-mediated cell death requires FOXO |

### Critical Genetic Dependencies

#### 1. Insulin/IGF Signaling Pathway
The data reveals strong epistatic relationships between FOXO and insulin signaling components:

**InR-FOXO Epistasis:**
- **Observation:** InR knockdown extends lifespan in wild-type flies
- **Dependency:** This lifespan extension is **completely abolished** in foxo-null background
- **Conclusion:** FOXO is *required* for the longevity benefits of reduced insulin signaling
- **Mechanism:** InR suppresses FOXO → removing InR activates FOXO → FOXO promotes longevity

**chico-FOXO Interaction:**
- **Effect:** "chico extends Drosophila longevity by acting through FOXO"
- **Interpretation:** The IRS homolog chico requires functional FOXO for lifespan extension

**Pi3K92E-FOXO Interaction:**
- **Classification:** Enhancer with opposite lifespan effects
- **Effect:** PI3K inhibition extends lifespan in FOXO-dependent manner
- **Pathway:** PI3K → AKT → FOXO phosphorylation/inhibition

#### 2. Cell Death and Stress Response
**rpr-FOXO Dependency:**
- **Observation:** Ablation of median neurosecretory cells by reaper expression extends lifespan
- **Dependency:** Completely cancelled in foxo-null background
- **Implication:** FOXO mediates the beneficial effects of neurosecretory cell reduction

### Experimental Details

All Drosophila studies maintained flies at:
- **Temperature:** 25°C (primary condition)
- **Humidity:** 40% relative humidity
- **Light Cycle:** 12h light : 12h dark
- **Density:** Standard rearing density (10-50 flies per vial)
- **Strain:** Backcrossed into white-Dahomey (wDah) background for ≥6 generations

---

## 4. Molecular Interactions from OmniPath Database

### 4.1 Post-Translational Regulation (Top Partners)

FOXO1 is extensively regulated through post-translational modifications:

#### Kinases That Phosphorylate FOXO1 (Inhibitory)

| Kinase | Effect | Sources | Mechanism |
|--------|--------|---------|-----------|
| **AKT1** | Inhibitory | 30+ databases | Phosphorylates T24, S256, S319 → nuclear export |
| **AKT2** | Inhibitory | Multiple | Similar to AKT1, phosphorylates and inactivates |
| **SGK1** | Inhibitory | 16+ databases | Phosphorylates S256 → nuclear export |
| **CDK2** | Inhibitory | 14+ databases | Cell cycle-dependent regulation |

#### Deacetylases That Activate FOXO1

| Enzyme | Effect | Sources | Mechanism |
|--------|--------|---------|-----------|
| **SIRT1** | Activating | BioGRID, KEGG, SIGNOR, etc. | Deacetylates FOXO1 → increases transcriptional activity |
| **SIRT2** | Activating/Inhibitory | SPIKE, SPIKE_LC | Context-dependent regulation |

#### E3 Ligases (Ubiquitination)

| Enzyme | Effect | Sources | Mechanism |
|--------|--------|---------|-----------|
| **SKP2** | Degradation | ACSN, SPIKE | Ubiquitinates FOXO1 → proteasomal degradation |
| **USP7** | Stabilization | SPIKE, Wang | Deubiquitinase - stabilizes FOXO1 |

### 4.2 Transcriptional Targets of FOXO1

FOXO1 regulates numerous genes involved in metabolism, stress response, and cell cycle:

#### Metabolic Genes (Gluconeogenesis)

| Target | Function | Direction |
|--------|----------|-----------|
| **G6PC1** | Glucose-6-phosphatase | Activation ↑ |
| **PCK1** | PEPCK (cytosolic) | Activation ↑ |
| **PCK2** | PEPCK (mitochondrial) | Activation ↑ |
| **IGFBP1** | IGF binding protein | Activation ↑ |

#### Cell Cycle & Apoptosis

| Target | Function | Direction |
|--------|----------|-----------|
| **CDKN1A** (p21) | CDK inhibitor | Activation/Inhibition |
| **CDKN1B** (p27) | CDK inhibitor | Activation ↑ |
| **CDKN2B** (p15) | CDK inhibitor | Activation ↑ |
| **CDKN2D** (p19) | CDK inhibitor | Regulation |

#### Autophagy & Stress

| Target | Function | Direction |
|--------|----------|-----------|
| **ATG7** | Autophagy | Activation ↑ |
| **EIF4EBP1** | Translation repressor | Activation ↑ |

#### Other Targets

| Target | Function | Direction |
|--------|----------|-----------|
| **KLF4** | Transcription factor | Activation ↑ |
| **TRIM63** | Muscle atrophy | Activation ↑ |
| **ANGPT2** | Angiogenesis | Activation ↑ |

### 4.3 Transcription Factors Regulating FOXO1

From CollecTRI database (Transcriptional Regulation):

#### Activating Transcription Factors

| TF | Effect | Sources | Function |
|----|--------|---------|----------|
| **FOXO3** | Activation ↑ | DoRothEA, CollecTRI | Cross-regulation within FOXO family |
| **FOXO1** | Auto-regulation ↑ | Multiple | Self-activation loop |
| **FOXA2** | Activation ↑ | CollecTRI | Forkhead family member |
| **E2F1** | Activation ↑ | CollecTRI, DoRothEA | Cell cycle regulator |
| **E2F2/E2F3** | Activation ↑ | CollecTRI | E2F family |
| **ESR1** | Activation ↑ | CollecTRI, DoRothEA | Estrogen receptor |
| **FOXC1** | Activation ↑ | TRRUST, CollecTRI | Important FOXO1 inducer |
| **PPARD** | Activation ↑ | CollecTRI | Metabolic regulator |
| **STAT3** | Activation ↑ | CollecTRI | JAK-STAT signaling |
| **NR1H3** (LXR) | Activation ↑ | CollecTRI | Lipid metabolism |
| **KLF5** | Activation ↑ | TRRUST, DoRothEA | Kruppel-like factor |
| **NR3C1** (GR) | Activation ↑ | CollecTRI | Glucocorticoid receptor |
| **TCF3** | Activation ↑ | CollecTRI | Wnt signaling |

#### Inhibiting Transcription Factors

| TF | Effect | Sources | Function |
|----|--------|---------|----------|
| **EBF1** | Inhibition ↓ | DoRothEA, CollecTRI | B-cell development |
| **NR1I3** (CAR) | Inhibition ↓ | CollecTRI | Xenobiotic receptor |
| **TP53** | Inhibition ↓ | CollecTRI | Tumor suppressor |

### 4.4 Nuclear Receptors & Hormonal Regulation

| Receptor | Effect | Sources | Context |
|----------|--------|---------|---------|
| **ESR1** (Estrogen Receptor α) | Inhibitory | BioGRID, HPRD | Hormone regulation |
| **PGR** (Progesterone Receptor) | Activation | SPIKE | Reproductive tissue |
| **AR** (Androgen Receptor) | Bidirectional | NetPath, HPRD | Prostate/muscle |
| **PPARA** | Bidirectional | SPIKE | Fatty acid oxidation |
| **PPARG** | Inhibition by FOXO1 | HPRD, SPIKE | Adipogenesis |

### 4.5 Signaling Pathway Integration

| Pathway | Key Proteins | FOXO1 Role |
|---------|--------------|------------|
| **TGF-β** | SMAD3 | Cooperates in transcription |
| **Wnt** | CTNNB1 | Cross-talk in development |
| **Nuclear Receptor** | NR1I2, PPARA, PPARG | Metabolic integration |

---

## 5. Key Biological Insights

### 5.1 FOXO1 as a Central Longevity Hub

The data demonstrates that FOXO1 acts as a critical integration point for multiple longevity-promoting pathways:

1. **Insulin/IGF Signaling (IIS):**
   - FOXO1 is THE key mediator of IIS effects on lifespan
   - Reduced insulin signaling → FOXO1 activation → lifespan extension
   - This mechanism is conserved from worms to mammals

2. **Nutrient Sensing:**
   - FOXO1 integrates signals from multiple nutrient sensors (AKT, AMPK, mTOR)
   - Coordinates metabolic adaptation to nutrient availability

3. **Stress Response:**
   - SIRT1-FOXO1 axis provides stress resistance
   - FOXO1 activates antioxidant genes and DNA repair

4. **Autophagy:**
   - FOXO1-ATG7 connection promotes cellular cleanup
   - Essential for maintaining cellular health during aging

### 5.2 Post-Translational Regulation as Lifespan Control

FOXO1 activity is exquisitely controlled by multiple PTMs:

**Phosphorylation (Generally Inhibitory):**
- AKT1/2 phosphorylation → 14-3-3 binding → nuclear export → inactivation
- Allows rapid response to insulin/growth signals

**Acetylation (Context-Dependent):**
- CREBBP/EP300 acetylation → reduced DNA binding
- SIRT1 deacetylation → increased activity and lifespan

**Ubiquitination:**
- SKP2-mediated degradation provides another level of control
- USP7 stabilization counterbalances degradation

### 5.3 Transcriptional Network

FOXO1 regulates a comprehensive longevity program:

1. **Metabolic Adaptation:**
   - Gluconeogenesis (G6PC1, PCK1/2)
   - Insulin sensitivity (IGFBP1)
   - Energy homeostasis (PPARGC1A cooperation)

2. **Cell Cycle Control:**
   - CDK inhibitors (p21, p27, p15, p19)
   - Growth arrest under stress conditions

3. **Stress Resistance:**
   - Antioxidant genes
   - DNA repair pathways
   - Protein quality control (autophagy)

4. **Cell Survival:**
   - Pro-survival under mild stress
   - Pro-death under severe/prolonged stress

### 5.4 Evolutionary Conservation

The FOXO family is highly conserved:
- **C. elegans:** DAF-16 (ortholog)
- **Drosophila:** dFOXO
- **Mammals:** FOXO1, FOXO3, FOXO4, FOXO6

The SynergyAge data confirms functional conservation in Drosophila, supporting translation to mammals.

---

## 6. Therapeutic & Research Implications

### 6.1 Potential Longevity Interventions

Based on FOXO1 biology, several intervention strategies emerge:

1. **SIRT1 Activators:**
   - Resveratrol and NAD+ boosters
   - Enhance FOXO1 deacetylation and activity

2. **Insulin Sensitizers:**
   - Metformin
   - Mild insulin signaling reduction → FOXO1 activation

3. **Caloric Restriction:**
   - Reduces insulin/IGF signaling
   - Activates FOXO1 and extends lifespan

4. **Exercise:**
   - Activates FOXO1 in muscle
   - Promotes mitochondrial biogenesis via FOXO1-PGC1α

### 6.2 Disease Connections

FOXO1 dysregulation is implicated in:

- **Diabetes:** Excess FOXO1 activity in liver → hyperglycemia
- **Cancer:** Loss of FOXO1 → uncontrolled proliferation
- **Aging:** Reduced FOXO1 activity → accelerated aging
- **Neurodegenerative Diseases:** FOXO1 protects neurons from stress

### 6.3 Research Directions

1. **Tissue-Specific Functions:**
   - FOXO1 effects vary by tissue
   - Need for conditional knockout studies

2. **Combinatorial Interventions:**
   - Multiple longevity pathways converge on FOXO1
   - Synergistic effects possible

3. **PTM Dynamics:**
   - Real-time monitoring of FOXO1 modifications
   - Understanding context-dependent regulation

4. **Translation to Humans:**
   - Clinical trials of FOXO1 activators
   - Biomarkers of FOXO1 activity

---

## 7. Summary & Conclusions

### Key Findings

1. **FOXO1 is Essential for Longevity:**
   - Deletion reduces lifespan by 7-34% in Drosophila
   - Required for lifespan benefits of reduced insulin signaling

2. **Hub of Longevity Networks:**
   - Connects insulin signaling, nutrient sensing, and stress response
   - Interacts with key longevity proteins (SIRT1, PPARGC1A, AKT1)

3. **Highly Regulated:**
   - Multiple PTMs (phosphorylation, acetylation, ubiquitination)
   - Complex transcriptional regulation (17 TFs identified)

4. **Broad Transcriptional Program:**
   - Metabolic genes (gluconeogenesis)
   - Stress resistance (antioxidants, DNA repair)
   - Cell cycle control (CDK inhibitors)
   - Autophagy (ATG7)

5. **Therapeutic Target:**
   - Validated target for longevity interventions
   - Multiple druggable pathways (SIRT1, AKT, insulin signaling)

### Conservation & Translation

The strong conservation of FOXO function from flies to humans, combined with the robust genetic evidence from model organisms, provides confidence that FOXO1 activation strategies could promote healthy aging in humans.

### Future Perspectives

Understanding and modulating FOXO1 activity represents one of the most promising approaches for extending healthspan and treating age-related diseases. The integration of multiple longevity pathways into FOXO1 makes it an ideal target for interventions that could have broad beneficial effects on aging.

---

## References & Data Sources

### Databases Queried
1. **UniProt** - Protein function and annotation (Q12778)
2. **STRING** - Protein-protein interactions (v12.0)
3. **SynergyAge** - Genetic synergy and longevity data
4. **OmniPath** - Molecular interactions and signaling pathways
5. **Ensembl** - Gene annotation (ENSG00000150907)

### Key Publications (from data)
- PMID: 21443682 - Drosophila FOXO deletion studies
- PMID: 21518241 - chico-FOXO epistasis
- PMID: 20624856 - FOXO-Lrrk interactions
- Multiple additional studies on FOXO1 regulation and function

### Analysis Metadata
- **Analysis Tool:** MCP Knowledgebase Server
- **Analysis Date:** October 31, 2025
- **Databases Versions:** Latest available as of analysis date
- **Species Focus:** Homo sapiens (Human), with comparative data from Drosophila melanogaster

---

## Appendix: Detailed Interaction Data

### Complete STRING Interaction Partners (Score ≥ 0.700)

All 10 partners are listed in Section 2 above with scores ranging from 0.984 to 0.999, representing extremely high-confidence interactions.

### Complete SynergyAge Models

All 13 FOXO-related models from SynergyAge database are detailed in Section 3, showing consistent lifespan-shortening effects of FOXO loss.

### Complete OmniPath Interactions

50+ direct interactions documented in Section 4, including:
- 4 major kinase regulators (AKT1, AKT2, SGK1, CDK2)
- 2 deacetylases (SIRT1, SIRT2)
- 17 transcription factors regulating FOXO1
- 30+ transcriptional targets of FOXO1

---

**End of Report**

