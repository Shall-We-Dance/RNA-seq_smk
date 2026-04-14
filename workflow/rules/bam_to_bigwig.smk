# rules/bam_to_bigwig.smk

OUTDIR = config["output"]["dir"]
KEEP_BAM = bool(config.get("output", {}).get("keep_bam", False))


def maybe_temp(path):
    return path if KEEP_BAM else temp(path)


rule bam_to_bigwig:
    """Generate raw and normalized bigWig coverage tracks from STAR BAM files."""
    input:
        bam = f"{OUTDIR}/star/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam"
    params:
        chrom_sizes = config.get("reference", {}).get(
            "chrom_sizes",
            config["reference"]["fasta"] + ".fai"
        ),
        bin_size = config.get("bigwig", {}).get("bin_size", 10),
        normalization = config.get("bigwig", {}).get("normalization", "RPKM"),
        exclude_flags = config.get("samtools_exclude_flags", 1804),
        remove_chrM_and_scaffolds = bool(
            config.get("bigwig", {}).get("remove_chrM_and_scaffolds", True)
        )
    output:
        raw_bw = f"{OUTDIR}/bigwig/{{sample}}/{{sample}}.raw.bw",
        norm_bw = f"{OUTDIR}/bigwig/{{sample}}/{{sample}}.normalized.bw",
        bam_index = maybe_temp(f"{OUTDIR}/star/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam.bai")
    log:
        f"logs/bigwig/{{sample}}.bamCoverage.log"
    threads: int(config["threads"]["deeptools"])
    conda:
        "envs/deeptools.yaml"
    shell:
        """
        GENOME_SIZE=$(awk '{{sum += $2}} END {{print sum}}' {params.chrom_sizes})
        mkdir -p $(dirname {output.raw_bw})

        # Ensure BAM index exists for bamCoverage.
        samtools index -@ {threads} "{input.bam}"

        EXTRA_BAMCOVERAGE_ARGS=""
        if [ "{params.remove_chrM_and_scaffolds}" = "True" ]; then
            EXCLUDE_BED=$(mktemp)
            awk 'BEGIN{{IGNORECASE=1}} {{
                chrom=$1; size=$2;
                if (chrom ~ /_/ || chrom ~ /scaffold|random|un|alt|fix|hap/ || chrom == "chrM" || chrom == "MT" || chrom == "M") {{
                    print chrom"\t0\t"size
                }}
            }}' {params.chrom_sizes} > "$EXCLUDE_BED"

            if [ -s "$EXCLUDE_BED" ]; then
                EXTRA_BAMCOVERAGE_ARGS="--blackListFileName $EXCLUDE_BED"
            fi
        fi
        
        # Unnormalized coverage
        bamCoverage \
            --bam {input.bam} \
            --outFileName {output.raw_bw} \
            --outFileFormat bigwig \
            --numberOfProcessors {threads} \
            --binSize {params.bin_size} \
            --normalizeUsing None \
            --samFlagExclude {params.exclude_flags} \
            $EXTRA_BAMCOVERAGE_ARGS \
            --effectiveGenomeSize $GENOME_SIZE \
            > {log} 2>&1
        
        # Normalized coverage
        bamCoverage \
            --bam {input.bam} \
            --outFileName {output.norm_bw} \
            --outFileFormat bigwig \
            --numberOfProcessors {threads} \
            --binSize {params.bin_size} \
            --normalizeUsing {params.normalization} \
            --samFlagExclude {params.exclude_flags} \
            $EXTRA_BAMCOVERAGE_ARGS \
            --effectiveGenomeSize $GENOME_SIZE \
            >> {log} 2>&1


        if [ -n "${{EXCLUDE_BED:-}}" ] && [ -f "$EXCLUDE_BED" ]; then
            rm -f "$EXCLUDE_BED"
        fi

        """
