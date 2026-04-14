# workflow/rules/align_star.smk
OUTDIR = config["output"]["dir"]


def maybe_temp(path):
    return path if KEEP_BAM else temp(path)



rule star_align:
    input:
        idx_sa=f"{config['reference']['star_index']}/SA",
        r1=f"{OUTDIR}/tmp/fastp/{{sample}}_R1.fastq.gz",
        r2=f"{OUTDIR}/tmp/fastp/{{sample}}_R2.fastq.gz"
    output:
        bam=maybe_temp(f"{OUTDIR}/star/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam"),
        log_final=f"{OUTDIR}/star/{{sample}}/{{sample}}.Log.final.out",
        log_final_qc=f"{OUTDIR}/qc/star/{{sample}}/{{sample}}.Log.final.out",
        sj=f"{OUTDIR}/star/{{sample}}/{{sample}}.SJ.out.tab"
    log:
        f"logs/star/{{sample}}.log"
    threads: int(config["threads"]["star"])
    conda:
        "envs/star.yaml"
    params:
        index=config["reference"]["star_index"],
        quant_mode="--quantMode GeneCounts" if config.get("star", {}).get("quant_mode_gene_counts", False) else ""
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.bam}) $(dirname {output.log_final}) $(dirname {output.log_final_qc}) $(dirname {log})

        STAR \
          --runThreadN {threads} \
          --genomeDir {params.index} \
          --readFilesIn {input.r1} {input.r2} \
          --readFilesCommand zcat \
          --outFileNamePrefix {OUTDIR}/star/{wildcards.sample}/{wildcards.sample}. \
          --outSAMtype BAM SortedByCoordinate \
          {params.quant_mode} \
          --outSJfilterOverhangMin 15 12 12 12 \
          --alignSJoverhangMin 15 \
          --alignSJDBoverhangMin 15 \
          --outFilterMultimapNmax 20 \
          --outFilterScoreMin 1 \
          --outFilterMatchNmin 1 \
          --outFilterMismatchNmax 2 \
          > {log} 2>&1

        cp {output.log_final} {output.log_final_qc}
        """
