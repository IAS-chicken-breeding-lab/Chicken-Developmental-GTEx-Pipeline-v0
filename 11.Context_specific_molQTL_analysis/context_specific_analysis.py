#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Context-Specific Expression and eGene Analysis Pipeline
-------------------------------------------------------
A modular, generalized pipeline to analyze context-specific (e.g., Sexes, Tissues,
Treatments) gene expression, eGene sharing, variant matching, and Linkage
Disequilibrium (LD) properties.

Author: GitHub Contributor
License: MIT
"""

import os
import glob
import itertools
import argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd


# ==============================================================================
# 1. Utility Functions
# ==============================================================================

def read_gene_list(file_path):
    """Read gene list from file (deduplicate and filter empty lines)."""
    genes = set()
    if os.path.exists(file_path):
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                gene = line.strip()
                if gene:
                    genes.add(gene)
    return genes


def write_gene_list(file_path, gene_set):
    """Write a set of genes to a specified file."""
    Path(os.path.dirname(file_path)).mkdir(parents=True, exist_ok=True)
    with open(file_path, 'w', encoding='utf-8') as f:
        for gene in sorted(gene_set):
            f.write(f"{gene}\n")


def parse_egene_file(file_path):
    """Parse eGene result file and extract pheno_id to variant_id mapping for eGene==yes."""
    gene_variant_map = {}
    if not os.path.exists(file_path):
        return gene_variant_map
    try:
        df = pd.read_csv(file_path, sep='\t')
        if 'pheno_id' in df.columns and 'variant_id' in df.columns:
            if 'eGene' in df.columns:
                df = df[df['eGene'] == 'yes']
            for _, row in df.iterrows():
                pheno_id = str(row['pheno_id']).strip()
                variant_id = str(row['variant_id']).strip()
                if pheno_id and variant_id and pd.notna(variant_id):
                    gene_variant_map[pheno_id] = variant_id
    except Exception as e:
        print(f"[Error] Failed to parse {file_path}: {e}")
    return gene_variant_map


# ==============================================================================
# 2. Pipeline Steps
# ==============================================================================

def step1_context_specific_expression_egene(input_dir, output_dir, contexts, groups):
    """Step 1: Context-specific and shared expression gene & eGene analysis."""
    for group in groups:
        grp_out_dir = os.path.join(output_dir, group) if group else output_dir
        dirs = {
            'spec': os.path.join(grp_out_dir, "1-1context-specific-expression-gene"),
            'share': os.path.join(grp_out_dir, "1-2context-share-expression-gene"),
            'spec_egene': os.path.join(grp_out_dir, "1-3context-specific-expression-egene")
        }
        for d in dirs.values():
            Path(d).mkdir(parents=True, exist_ok=True)

        ctx_genes = {}
        ctx_egenes = {}
        for c in contexts:
            fname = f"{c}_{group}.txt" if group else f"{c}.txt"
            egene_fname = f"{c}_{group}_eGene.txt" if group else f"{c}_eGene.txt"
            
            ctx_genes[c] = read_gene_list(os.path.join(input_dir, "expression-gene", fname))
            ctx_egenes[c] = read_gene_list(os.path.join(input_dir, "eGene", egene_fname))

        ctx_stats = {c: {'total': len(ctx_egenes[c]), 'union': set(), 'list': []} for c in contexts}
        all_egene_sum = sum(len(v) for v in ctx_egenes.values())

        for c1, c2 in itertools.combinations(contexts, 2):
            g1, g2 = ctx_genes[c1], ctx_genes[c2]
            spec_1, spec_2 = g1 - g2, g2 - g1
            shared = g1 & g2

            write_gene_list(os.path.join(dirs['spec'], f"{c1}-{c2}.txt"), spec_1)
            write_gene_list(os.path.join(dirs['spec'], f"{c2}-{c1}.txt"), spec_2)
            write_gene_list(os.path.join(dirs['share'], f"{c1}-{c2}.txt"), shared)
            write_gene_list(os.path.join(dirs['share'], f"{c2}-{c1}.txt"), shared)

            spec_eg1, spec_eg2 = spec_1 & ctx_egenes[c1], spec_2 & ctx_egenes[c2]
            write_gene_list(os.path.join(dirs['spec_egene'], f"{c1}-{c2}.txt"), spec_eg1)
            write_gene_list(os.path.join(dirs['spec_egene'], f"{c2}-{c1}.txt"), spec_eg2)

            for c_curr, eg_set in [(c1, spec_eg1), (c2, spec_eg2)]:
                ctx_stats[c_curr]['union'].update(eg_set)
                ctx_stats[c_curr]['list'].extend(list(eg_set))

        summary_file = os.path.join(grp_out_dir, "1_summary_statistics.tsv")
        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write("Context\tTotal_eGene\tUnique_Specific_eGene\tPct_In_Context(%)\tTotal_Specific_eGene\tPct_In_All(%)\n")
            for c in contexts:
                u_cnt = len(ctx_stats[c]['union'])
                u_pct = (u_cnt / ctx_stats[c]['total'] * 100) if ctx_stats[c]['total'] > 0 else 0
                l_cnt = len(ctx_stats[c]['list'])
                l_pct = (l_cnt / all_egene_sum * 100) if all_egene_sum > 0 else 0
                f.write(f"{c}\t{ctx_stats[c]['total']}\t{u_cnt}\t{u_pct:.2f}%\t{l_cnt}\t{l_pct:.2f}%\n")


def step2_context_specific_egene_share_gene(input_dir, output_dir, contexts, groups):
    """Step 2: Intersection analysis of eGenes and shared expressed genes."""
    for group in groups:
        grp_out_dir = os.path.join(output_dir, group) if group else output_dir
        dirs = {
            'spec_egene': os.path.join(grp_out_dir, "2-1context-specific-egene"),
            'share_egene': os.path.join(grp_out_dir, "2-2context-share-egene"),
            'spec_egene_share_gene': os.path.join(grp_out_dir, "2-3context-specific-egene_share-gene"),
            'share_egene_share_gene': os.path.join(grp_out_dir, "2-4context-share-egene_share-gene")
        }
        for d in dirs.values():
            Path(d).mkdir(parents=True, exist_ok=True)

        ctx_egenes = {}
        for c in contexts:
            egene_fname = f"{c}_{group}_eGene.txt" if group else f"{c}_eGene.txt"
            ctx_egenes[c] = read_gene_list(os.path.join(input_dir, "eGene", egene_fname))

        for c1, c2 in itertools.combinations(contexts, 2):
            eg1, eg2 = ctx_egenes[c1], ctx_egenes[c2]
            spec_eg1, spec_eg2 = eg1 - eg2, eg2 - eg1
            share_eg = eg1 & eg2

            write_gene_list(os.path.join(dirs['spec_egene'], f"{c1}-{c2}.txt"), spec_eg1)
            write_gene_list(os.path.join(dirs['spec_egene'], f"{c2}-{c1}.txt"), spec_eg2)
            write_gene_list(os.path.join(dirs['share_egene'], f"{c1}-{c2}.txt"), share_eg)
            write_gene_list(os.path.join(dirs['share_egene'], f"{c2}-{c1}.txt"), share_eg)

            share_expr_file = os.path.join(grp_out_dir, "1-2context-share-expression-gene", f"{c1}-{c2}.txt")
            if os.path.exists(share_expr_file):
                share_expr_set = read_gene_list(share_expr_file)
                write_gene_list(os.path.join(dirs['spec_egene_share_gene'], f"{c1}-{c2}.txt"), spec_eg1 & share_expr_set)
                write_gene_list(os.path.join(dirs['spec_egene_share_gene'], f"{c2}-{c1}.txt"), spec_eg2 & share_expr_set)
                
                share_eg_share_gene = share_eg & share_expr_set
                write_gene_list(os.path.join(dirs['share_egene_share_gene'], f"{c1}-{c2}.txt"), share_eg_share_gene)
                write_gene_list(os.path.join(dirs['share_egene_share_gene'], f"{c2}-{c1}.txt"), share_eg_share_gene)


def step3_1_variant_match(input_dir, output_dir, egene_results_dir, contexts, groups):
    """Step 3-1: Extract Variant IDs for shared eGenes."""
    for group in groups:
        grp_out_dir = os.path.join(output_dir, group) if group else output_dir
        share_gene_dir = os.path.join(grp_out_dir, "2-4context-share-egene_share-gene")
        out_dir = os.path.join(grp_out_dir, "3-1context-share-egene_share-gene-variant")

        if not os.path.exists(share_gene_dir):
            continue
        Path(out_dir).mkdir(parents=True, exist_ok=True)

        ctx_variant_maps = {}
        for c in contexts:
            egene_file = os.path.join(egene_results_dir, f"{c}_{group}_LMM_eGene.txt" if group else f"{c}_LMM_eGene.txt")
            ctx_variant_maps[c] = parse_egene_file(egene_file)

        for fname in os.listdir(share_gene_dir):
            if fname.endswith('.txt'):
                basename = fname.replace('.txt', '')
                if '-' not in basename:
                    continue
                c1, c2 = basename.split('-')[:2]

                share_genes = read_gene_list(os.path.join(share_gene_dir, fname))
                results = []
                for gene in share_genes:
                    v1 = ctx_variant_maps.get(c1, {}).get(gene, "NA")
                    v2 = ctx_variant_maps.get(c2, {}).get(gene, "NA")
                    if v1 != "NA" or v2 != "NA":
                        results.append([gene, v1, v2])

                if results:
                    df = pd.DataFrame(results, columns=['pheno_id', f'variant_id_{c1}', f'variant_id_{c2}'])
                    df.to_csv(os.path.join(out_dir, fname), sep='\t', index=False)


def step3_2_variant_merge(output_dir, groups):
    """Step 3-2: Merge and standardize deduplication of variant pairs."""
    for group in groups:
        grp_out_dir = os.path.join(output_dir, group) if group else output_dir
        in_dir = os.path.join(grp_out_dir, "3-1context-share-egene_share-gene-variant")
        out_dir = os.path.join(grp_out_dir, "3-2context-share-egene_share-gene-variant-LD")
        out_file = os.path.join(out_dir, "all-contexts-eGene-variant-LD.txt")

        if not os.path.exists(in_dir):
            continue
        Path(out_dir).mkdir(parents=True, exist_ok=True)

        dfs = []
        for fname in os.listdir(in_dir):
            if fname.endswith('.txt') and not fname.startswith('.'):
                try:
                    df = pd.read_csv(os.path.join(in_dir, fname), sep='\t')
                    if len(df.columns) >= 3 and df.columns[0] == 'pheno_id':
                        df = df.rename(columns={df.columns[1]: 'variant_id_1', df.columns[2]: 'variant_id_2'})
                        dfs.append(df[['pheno_id', 'variant_id_1', 'variant_id_2']])
                except Exception:
                    pass

        if dfs:
            combined = pd.concat(dfs, ignore_index=True).fillna('NA')
            sorted_vars = combined.apply(lambda r: pd.Series(sorted([str(r['variant_id_1']), str(r['variant_id_2'])])), axis=1)
            sorted_vars.columns = ['v1_s', 'v2_s']
            
            temp = pd.concat([combined, sorted_vars], axis=1)
            dedup = temp.drop_duplicates(subset=['pheno_id', 'v1_s', 'v2_s'])
            final_df = dedup[['pheno_id', 'variant_id_1', 'variant_id_2']].sort_values('pheno_id')
            final_df.to_csv(out_file, sep='\t', index=False)


def _process_ld_worker(args):
    """Worker logic for parallel LD matching on a single file."""
    file_path, output_dir, ld_dict = args
    fname = Path(file_path).name
    try:
        lines = []
        with open(file_path, 'r', encoding='utf-8') as f:
            for i, line in enumerate(f):
                parts = line.strip().split('\t')
                if i == 0 and parts[0] == 'pheno_id':
                    lines.append('\t'.join(parts + ['LD']))
                    continue
                if len(parts) >= 3:
                    p_id, v1, v2 = parts[0].strip(), parts[1].strip(), parts[2].strip()
                    ld_val = ld_dict.get(f"{p_id}_{v1}_{v2}", ld_dict.get(f"{p_id}_{v2}_{v1}", 'NA'))
                    lines.append('\t'.join(parts + [ld_val]))
        
        with open(os.path.join(output_dir, f"matched_{fname}"), 'w', encoding='utf-8') as f:
            f.write('\n'.join(lines) + '\n')
    except Exception as e:
        print(f"[Error] Failed matching {fname}: {e}")


def step3_3_ld_match(output_dir, groups, max_workers=4):
    """Step 3-3: Parallel mapping of calculated LD results."""
    for group in groups:
        grp_out_dir = os.path.join(output_dir, group) if group else output_dir
        in_dir = os.path.join(grp_out_dir, "3-1context-share-egene_share-gene-variant")
        ld_res_file = os.path.join(grp_out_dir, "3-2context-share-egene_share-gene-variant-LD", "all-contexts-eGene-variant-LD_ld_result.txt")
        out_dir = os.path.join(grp_out_dir, "3-2context-share-egene_share-gene-variant-LD")

        if not os.path.exists(in_dir) or not os.path.exists(ld_res_file):
            continue

        ld_df = pd.read_csv(ld_res_file, sep='\t', header=None, dtype=str)
        ld_dict = {}
        for row in ld_df.values:
            if len(row) >= 4:
                p_id, v1, v2, ld = str(row[0]).strip(), str(row[1]).strip(), str(row[2]).strip(), str(row[3]).strip()
                ld_dict[f"{p_id}_{v1}_{v2}"] = ld
                ld_dict[f"{p_id}_{v2}_{v1}"] = ld

        file_list = glob.glob(os.path.join(in_dir, "*.txt"))
        args_list = [(fpath, out_dir, ld_dict) for fpath in file_list]

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(_process_ld_worker, arg) for arg in args_list]
            for future in as_completed(futures):
                future.result()


def step3_4_ld_filter_and_summary(input_dir, output_dir, contexts, groups):
    """Step 3-4: Filter LD thresholds and generate global summary report."""
    summary_rows = []

    for group in groups:
        grp_out_dir = os.path.join(output_dir, group) if group else output_dir
        matched_dir = os.path.join(grp_out_dir, "3-2context-share-egene_share-gene-variant-LD")
        dir_ld02 = os.path.join(grp_out_dir, "3-3context-share-egene_share-gene-variant-LD0.2")
        dir_ld08 = os.path.join(grp_out_dir, "3-4context-share-egene_share-gene-variant-LD0.8")

        Path(dir_ld02).mkdir(parents=True, exist_ok=True)
        Path(dir_ld08).mkdir(parents=True, exist_ok=True)

        for fpath in glob.glob(os.path.join(matched_dir, "matched_*.txt")):
            fname = os.path.basename(fpath)
            try:
                df = pd.read_csv(fpath, sep='\t')
                if 'LD' in df.columns:
                    df['LD_num'] = pd.to_numeric(df['LD'], errors='coerce')
                    df[df['LD_num'] < 0.2].drop(columns=['LD_num']).to_csv(os.path.join(dir_ld02, fname), sep='\t', index=False)
                    df[df['LD_num'] > 0.8].drop(columns=['LD_num']).to_csv(os.path.join(dir_ld08, fname), sep='\t', index=False)
            except Exception:
                pass

        for ctx in contexts:
            egene_fname = f"{ctx}_{group}_eGene.txt" if group else f"{ctx}_eGene.txt"
            total_egenes = len(read_gene_list(os.path.join(input_dir, "eGene", egene_fname)))

            genes_ld02, genes_ld08 = set(), set()
            for fpath in glob.glob(os.path.join(dir_ld02, "matched_*.txt")):
                if f"-{ctx}.txt" in fpath or f"matched_{ctx}-" in fpath:
                    genes_ld02.update(read_gene_list(fpath))

            for fpath in glob.glob(os.path.join(dir_ld08, "matched_*.txt")):
                if f"-{ctx}.txt" in fpath or f"matched_{ctx}-" in fpath:
                    genes_ld08.update(read_gene_list(fpath))

            c_ld02, c_ld08 = len(genes_ld02), len(genes_ld08)
            summary_rows.append({
                "Group": group if group else "All",
                "Context": ctx,
                "Total_eGene": total_egenes,
                "LD_lt_0.2_Count": c_ld02,
                "LD_lt_0.2_Pct": f"{(c_ld02/total_egenes*100):.2f}%" if total_egenes > 0 else "0.00%",
                "LD_gt_0.8_Count": c_ld08,
                "LD_gt_0.8_Pct": f"{(c_ld08/total_egenes*100):.2f}%" if total_egenes > 0 else "0.00%"
            })

    if summary_rows:
        summary_df = pd.DataFrame(summary_rows)
        out_summary = os.path.join(output_dir, "ld_extract_summary.tsv")
        summary_df.to_csv(out_summary, sep='\t', index=False)


# ==============================================================================
# 3. Main Entrypoint
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Context-Specific Gene Expression & eGene Analysis Pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("-i", "--input_dir", type=str, required=True, help="Input directory containing expression and eGene lists")
    parser.add_argument("-o", "--output_dir", type=str, required=True, help="Output directory for results")
    parser.add_argument("-e", "--egene_results_dir", type=str, default=None, help="Directory containing raw eQTL/LMM results (default: <input_dir>/eGene_results)")
    parser.add_argument("-c", "--contexts", nargs="+", required=True, help="List of contexts to compare (e.g. 'Male Female' or 'TissueA TissueB')")
    parser.add_argument("-g", "--groups", nargs="*", default=[""], help="Optional subgroup tags (e.g., tissue IDs if comparing sexes per tissue)")
    parser.add_argument("-s", "--steps", nargs="+", default=['all'], choices=['1', '2', '3-1', '3-2', '3-3', '3-4', 'all'], help="Pipeline steps to execute")
    parser.add_argument("-w", "--workers", type=int, default=4, help="Parallel worker threads for LD matching")

    args = parser.parse_args()

    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output_dir)
    egene_results_dir = os.path.abspath(args.egene_results_dir) if args.egene_results_dir else os.path.join(input_dir, "eGene_results")

    steps = set(args.steps)
    if 'all' in steps:
        steps = {'1', '2', '3-1', '3-2', '3-3', '3-4'}

    groups = args.groups if args.groups else [""]

    if '1' in steps:
        step1_context_specific_expression_egene(input_dir, output_dir, args.contexts, groups)
    if '2' in steps:
        step2_context_specific_egene_share_gene(input_dir, output_dir, args.contexts, groups)
    if '3-1' in steps:
        step3_1_variant_match(input_dir, output_dir, egene_results_dir, args.contexts, groups)
    if '3-2' in steps:
        step3_2_variant_merge(output_dir, groups)
    if '3-3' in steps:
        step3_3_ld_match(output_dir, groups, max_workers=args.workers)
    if '3-4' in steps:
        step3_4_ld_filter_and_summary(input_dir, output_dir, args.contexts, groups)

    print("Pipeline processing completed successfully.")


if __name__ == "__main__":
    main()