import os
import pandas as pd
import requests
import time
from typing import Dict, Any, List, Optional, Tuple
from dotenv import load_dotenv
from alphagenome.data import genome
from alphagenome.models import dna_client
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

load_dotenv()

ALPHA_GENOME_API_KEY = os.getenv("ALPHA_GENOME_API_KEY")
ENSEMBL_REST_URL = "https://rest.ensembl.org"

# Chromosome mapping from NCBI accession to UCSC format
CHROMOSOME_MAPPING = {
    'NC_000001.11': 'chr1', 'NC_000002.12': 'chr2', 'NC_000003.12': 'chr3',
    'NC_000004.12': 'chr4', 'NC_000005.10': 'chr5', 'NC_000006.12': 'chr6',
    'NC_000007.14': 'chr7', 'NC_000008.11': 'chr8', 'NC_000009.12': 'chr9',
    'NC_000010.11': 'chr10', 'NC_000011.10': 'chr11', 'NC_000012.12': 'chr12',
    'NC_000013.11': 'chr13', 'NC_000014.9': 'chr14', 'NC_000015.10': 'chr15',
    'NC_000016.10': 'chr16', 'NC_000017.11': 'chr17', 'NC_000018.10': 'chr18',
    'NC_000019.10': 'chr19', 'NC_000020.11': 'chr20', 'NC_000021.9': 'chr21',
    'NC_000022.11': 'chr22', 'NC_000023.11': 'chrX', 'NC_000024.10': 'chrY',
}

def get_ensembl_sequence(chrom: str, start: int, end: int, strand: int = 1) -> Optional[str]:
    """Fetch reference sequence from Ensembl REST API (GRCh38)."""
    chrom_clean = CHROMOSOME_MAPPING.get(chrom, chrom)
    if chrom_clean.startswith('chr'):
        chrom_clean = chrom_clean[3:]
    elif chrom_clean.startswith('NC_'):
        for nc, chr_name in CHROMOSOME_MAPPING.items():
            if nc == chrom:
                chrom_clean = chr_name[3:]
                break
    
    fetch_start = min(start, end)
    fetch_end = max(start, end)

    url = f"{ENSEMBL_REST_URL}/sequence/region/human/{chrom_clean}:{fetch_start}..{fetch_end}:{strand}?content-type=application/json"
    try:
        response = requests.get(url, headers={"Content-Type": "application/json"})
        response.raise_for_status()
        data = response.json()
        return data.get("seq")
    except Exception as e:
        logger.error(f"Error fetching sequence for {chrom}:{start}-{end}: {e}")
        return None

def get_snp_info_ensembl(rsid: str) -> Optional[Dict[str, Any]]:
    """Fetch detailed SNP info from Ensembl Variation API (GRCh38)."""
    url = f"{ENSEMBL_REST_URL}/variation/human/{rsid}?content-type=application/json"
    try:
        response = requests.get(url, headers={"Content-Type": "application/json"})
        response.raise_for_status()
        data = response.json()
        
        mappings = data.get('mappings', [])
        for mapping in mappings:
            if mapping.get('assembly_name') == 'GRCh38':
                allele_string = mapping.get('allele_string', '')
                if '/' in allele_string:
                    ref = allele_string.split('/')[0]
                    alt = allele_string.split('/')[1]
                    return {
                        'chrom': f"chr{mapping['seq_region_name']}",
                        'pos': mapping['start'],
                        'ref': ref,
                        'alt': alt,
                        'strand': mapping['strand']
                    }
        return None
    except Exception as e:
        return None

def save_fasta(filename: str, header: str, sequence: str):
    """Save a sequence to a FASTA file."""
    with open(filename, 'w') as f:
        f.write(f">{header}\n")
        for i in range(0, len(sequence), 60):
            f.write(sequence[i:i+60] + "\n")

def run_sequence_comparer():
    logger.info("Reading input data from overlapping_genes_with_snps.xlsx")
    df = pd.read_excel('overlapping_genes_with_snps.xlsx')
    
    # Prioritize significant SNPs, then others until we reach 70
    sig_df = df[df['SNP Association'] == 'significant'].copy()
    non_sig_df = df[df['SNP Association'] != 'significant'].copy()
    
    # We'll iterate through sig_df first, then non_sig_df
    all_potential_rows = pd.concat([sig_df, non_sig_df])
    
    logger.info(f"Scanning through {len(all_potential_rows)} potential Gene/SNP pairs to find 70 valid ones.")

    os.makedirs('data/fasta', exist_ok=True)

    model = None
    if ALPHA_GENOME_API_KEY:
        try:
            model = dna_client.create(ALPHA_GENOME_API_KEY)
            logger.info("AlphaGenome client initialized.")
        except Exception as e:
            logger.error(f"Failed to initialize AlphaGenome client: {e}")

    gene_sequences = {}
    results = []
    processed_count = 0
    window = 5000 

    for idx, row in all_potential_rows.iterrows():
        gene_symbol = row['Gene Symbol']
        snp_id = row['SNP Identifier']
        
        if not str(snp_id).startswith('rs'):
            continue

        snp_info = get_snp_info_ensembl(snp_id)
        if not snp_info:
            continue

        if gene_symbol not in gene_sequences:
            chrom = row['Chromosome']
            start = int(row['Start'])
            end = int(row['End'])
            
            w_start = min(start, end) - window
            w_end = max(start, end) + window
            
            ref_seq = get_ensembl_sequence(chrom, w_start, w_end)
            if ref_seq:
                gene_sequences[gene_symbol] = {
                    'seq': ref_seq,
                    'chrom': chrom,
                    'start': w_start,
                    'end': w_end
                }
                save_fasta(f"data/fasta/{gene_symbol}_ref.fasta", f"{gene_symbol}_ref", ref_seq)
            else:
                continue

        gene_info = gene_sequences[gene_symbol]
        snp_pos = snp_info['pos']
        rel_pos = snp_pos - gene_info['start']
        
        if rel_pos < 0 or rel_pos >= len(gene_info['seq']):
            continue

        mut_seq_list = list(gene_info['seq'])
        # Handle cases where alt might be longer than 1bp (indels) - AlphaGenome usually works with SNPs
        if len(snp_info['ref']) == 1 and len(snp_info['alt']) == 1:
            mut_seq_list[rel_pos] = snp_info['alt']
            mut_seq = "".join(mut_seq_list)
            
            alt_fasta_path = f"data/fasta/{gene_symbol}_{snp_id}_alt.fasta"
            save_fasta(alt_fasta_path, f"{gene_symbol}_{snp_id}_alt", mut_seq)

            if model:
                try:
                    # Use the closest supported sequence length if possible, or truncate/extend
                    # Supported lengths: [2048, 16384, 131072, 524288, 1048576]
                    # We will use 131072 (128kb) as a default if the gene is smaller, 
                    # or 1048576 (1MB) if it's larger but within bounds.
                    # However, for the 'interval' argument, we should probably just center 
                    # a window of fixed size around the SNP.
                    
                    target_length = 131072 # 128kb
                    half_len = target_length // 2
                    
                    # Define a fixed-size interval centered on the SNP
                    ag_chrom = CHROMOSOME_MAPPING.get(gene_info['chrom'], gene_info['chrom'])
                    if not ag_chrom.startswith('chr'): ag_chrom = f"chr{ag_chrom}"
                    
                    fixed_start = snp_pos - half_len
                    fixed_end = snp_pos + half_len
                    
                    interval = genome.Interval(ag_chrom, fixed_start, fixed_end)
                    variant = genome.Variant(ag_chrom, snp_pos, snp_info['ref'], snp_info['alt'])
                    
                    outputs = model.predict_variant(
                        interval=interval,
                        variant=variant,
                        ontology_terms=['UBERON:0001157'], # Caudate nucleus (Brain)
                        requested_outputs=[dna_client.OutputType.RNA_SEQ]
                    )
                    results.append({'gene': gene_symbol, 'snp': snp_id, 'status': 'success', 'outputs': outputs})
                except Exception as e:
                    results.append({'gene': gene_symbol, 'snp': snp_id, 'status': 'failed', 'error': str(e)})
            else:
                results.append({'gene': gene_symbol, 'snp': snp_id, 'status': 'skipped', 'reason': 'No API Key'})
            
            processed_count += 1
            logger.info(f"Successfully processed {processed_count}/70: {gene_symbol} {snp_id}")
            
            if processed_count >= 70:
                break
        
        time.sleep(0.05)

    results_df = pd.DataFrame(results)
    results_df.to_csv('alphagenome_comparison_results.csv', index=False)
    logger.info(f"Analysis complete. Processed {processed_count} Gene/SNP pairs. Results saved to alphagenome_comparison_results.csv")

if __name__ == "__main__":
    run_sequence_comparer()
