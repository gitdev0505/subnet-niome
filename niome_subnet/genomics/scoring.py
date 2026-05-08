#!/usr/bin/env python3
import bittensor as bt
import pysam
import subprocess
import os

from collections import defaultdict
from niome_subnet.genomics.model import GroundTruth, MinerScore, MinerSubmission
from niome_subnet.utils.constants import FORWARD_TIMEOUT

# -----------------------------
# Compress and index vcf
# -----------------------------
def preprocess_vcf(vcf_path: str) -> str:
    subprocess.run(f"bgzip -c {vcf_path} > {vcf_path}.gz", shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    subprocess.run(f"tabix -f -p vcf {vcf_path}.gz", shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return f"{vcf_path}.gz"


# -----------------------------
# Create mapping file
# -----------------------------
def create_mapping_file(ref_fasta: str, read1: str, read2: str) -> str:
    bam_path = "data/sim.bam"
    if os.path.exists(bam_path):
        return bam_path

    # bwa index
    subprocess.run(f"bwa index {ref_fasta}", shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # align
    subprocess.run(f"bwa mem {ref_fasta} {read1} {read2} | samtools sort -o {bam_path}", shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    subprocess.run(f"samtools index {bam_path}", shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    return bam_path


# -----------------------------
# 1. VCf LOADER
# -----------------------------
def load_vcf(path):
    vcf = pysam.VariantFile(path)
    variants = set()

    for rec in vcf.fetch():
        if rec.alts is None:
            continue
        for alt in rec.alts:
            variants.add((rec.contig, rec.pos, rec.ref, alt))

    return variants


# -----------------------------
# 2. NORMALIZATION VIA BCFTOOLS
# -----------------------------
def normalize_vcf(vcf_in, ref_fai, out):
    subprocess.run([
        "bcftools", "norm",
        "-f", ref_fai,
        "-c", "x",
        "-m", "-both",
        vcf_in,
        "-Oz",
        "-o", out
    ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    subprocess.run(["bcftools", "index", "-f", out], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return out


# -----------------------------
# 3. LOAD DEPTH FROM BAM
# -----------------------------
def load_depth(bam_path):
    bam = pysam.AlignmentFile(bam_path, "rb")
    depth = defaultdict(int)

    for pileup_col in bam.pileup():
        depth[pileup_col.reference_pos + 1] = pileup_col.nsegments

    return depth


# -----------------------------
# 4. COMPUTE CONFUSION MATRIX
# -----------------------------
def compute_sets(truth, pred):
    tp = truth & pred
    fp = pred - truth
    fn = truth - pred
    return tp, fp, fn


# -----------------------------
# 5. VARIANT DIFFICULTY MODEL
# -----------------------------
def variant_weight(depth):
    # Simple but effective heuristic
    if depth < 10:
        return 0.3   # hard
    elif depth < 20:
        return 0.6   # medium
    else:
        return 1.0   # easy


# -----------------------------
# 6. WEIGHTED METRICS
# -----------------------------
def weighted_metrics(tp, fp, fn, depth_map):
    tp_w = sum(variant_weight(depth_map.get(v[1], 0)) for v in tp)
    fp_w = sum(variant_weight(depth_map.get(v[1], 0)) for v in fp)
    fn_w = sum(variant_weight(depth_map.get(v[1], 0)) for v in fn)

    precision = tp_w / (tp_w + fp_w + 1e-9)
    recall = tp_w / (tp_w + fn_w + 1e-9)
    f1 = 2 * precision * recall / (precision + recall + 1e-9)

    return precision, recall, f1


# -----------------------------
# 7. FINAL SCORE
# -----------------------------
def final_score(p, r, f1, response_time):
    score1 = 0.4 * f1 + 0.3 * p + 0.3 * r
    # score2 = (max(FORWARD_TIMEOUT - response_time, 0) / FORWARD_TIMEOUT) ** 2 * score1
    return score1


# -----------------------------
# 8. MAIN PIPELINE
# -----------------------------
def score(miner_submission: MinerSubmission, ground_truth: GroundTruth, bam: str) -> MinerScore:
    try:
        miner_origin_vcf = "data/miner.vcf"
        with open(miner_origin_vcf, "w") as f:
            f.write(miner_submission.vcf_content)

        miner_vcf = preprocess_vcf(miner_origin_vcf)

        # [1] Normalizing VCFs
        truth_norm = "data/truth.norm.vcf.gz"
        miner_norm = "data/miner.norm.vcf.gz"

        normalize_vcf(ground_truth.truth_vcf, ground_truth.ref, truth_norm)
        normalize_vcf(miner_vcf, ground_truth.ref, miner_norm)

        # [2] Loading variants
        truth = load_vcf(truth_norm)
        pred  = load_vcf(miner_norm)

        # [3] Loading depth from BAM
        depth = load_depth(bam)

        # [4] Computing overlap
        tp, fp, fn = compute_sets(truth, pred)

        # [5] Computing weighted metrics
        p, r, f1 = weighted_metrics(tp, fp, fn, depth)

        score_val = final_score(p, r, f1, miner_submission.response_time)

        miner_score = MinerScore(
            uid=miner_submission.uid,
            precision=p,
            recall=r,
            f1_score=f1,
            response_time=miner_submission.response_time,
            final_score=score_val,
            log=f"Scored with {len(tp)} TP, {len(fp)} FP, {len(fn)} FN.\nPrecision: {p:.4f}, Recall: {r:.4f}, F1: {f1:.4f}, Final Score: {score_val:.4f}\n\nMiner VCF\n{miner_submission.vcf_content}",
        )

        return miner_score
    except Exception as e:
        bt.logging.warning(f"Error scoring miner {miner_submission.uid}: {e}")
        return MinerScore(
            uid=miner_submission.uid,
            precision=0.0,
            recall=0.0,
            f1_score=0.0,
            response_time=miner_submission.response_time,
            final_score=0.0,
            log=f"Error: {e}",
        )
