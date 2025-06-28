Single#!/bin/bash
# filepath: /home/fede118/Arabidopsis_MIRAST_WGCNA/run_wgcna_pipeline.sh

# Darle permisos de ejecución: 
# chmod +x run_wgcna_pipeline.sh

# Como ejecutar:
# ./run_wgcna_pipeline.sh

# Deben estar instalados previamente los programas fastq-dump (v3.0.2) de SRAtoolkit, bowtie(v1.3.1), hisat2(v2.2.1), fastQC (v0.12.1), samtools (v1.20), stringtie(v2.2.1).

# Configuración de rutas. El orden de los directorios se encuentra en el mismo directorio que este script. "WORKDIR==/Arabidopsis_MIRAST_WGCNA"
WORKDIR=$(pwd)
FASTQ_DIR="$WORKDIR/fastq"
FILTERED_DIR="$WORKDIR/filtered_reads"
FASTQC_DIR="$WORKDIR/fastQC_reports"
HISAT_OUT="$WORKDIR/hisat_out"
STRINGTIE_OUT="$WORKDIR/stringtie_out"
GENOMIC_IDX="$WORKDIR/genomic_files/TAIR10/TAIR10idx"
GFF="$WORKDIR/genomic_files/Arapor11/Araport11_GFF3_genes_transposons.20241001_modifed_pri-MIR394b.gff"
RRNA_IDX="$WORKDIR/rRNA_Arabidopsis/At_rRNA"

# Crear directorios si no existen
mkdir -p "$FASTQ_DIR" "$FILTERED_DIR/reports" "$FASTQC_DIR" "$HISAT_OUT/reports" "$STRINGTIE_OUT/gene_abundance"

# Los directorios rRNA_Arabidopsis y genomic_files deben tener los archivos de rRNA dos subunidades (con U reemplazado por T en las secuencias) indexados y los directorios Arapor11 (con el archivo .gtf utilizado) y TAIR10 con el genoma de Arabidopsis indexado, respectivamente.

# Leer lista de muestras: chequeo de prescencia del archivo con la lista de los códigos de Run de SRA
# Formato: SRRxxx Single|Paired
SAMPLES_FILE="samples.txt"
if [ ! -f "$SAMPLES_FILE" ]; then
    echo "Crea un archivo samples.txt con una muestra por línea: SRRxxx Single|Paired NA|F|R|RF"
    exit 1
fi

while read -r SRR TYPE STRAND; do
    echo "Procesando $SRR ($TYPE-end)"

    # 1. Descargar fastq
    if [ "$TYPE" == "Single" ]; then
        fastq-dump --gzip --outdir "$FASTQ_DIR" "$SRR"
        FASTQ="$FASTQ_DIR/${SRR}.fastq.gz"
    else
        fastq-dump --split-3 --gzip --outdir "$FASTQ_DIR" "$SRR"
        FASTQ1="$FASTQ_DIR/${SRR}_1.fastq.gz"
        FASTQ2="$FASTQ_DIR/${SRR}_2.fastq.gz"
    fi

    # 2. Control de calidad con FastQC
    if [ "$TYPE" == "Single" ]; then
        fastqc -o "$FASTQC_DIR" "$FASTQ"
    else
        fastqc -o "$FASTQC_DIR" "$FASTQ1" "$FASTQ2"
    fi

    echo "Revisa los reportes de calidad en $FASTQC_DIR y presiona [Enter] para continuar o Ctrl+C para abortar."
    read

    # 3. Filtrado de rRNA con bowtie
    if [ "$TYPE" == "Single" ]; then
        bowtie --quiet -v 3 -x "$RRNA_IDX" -q "$FASTQ" --un "$FILTERED_DIR/${SRR}_clean.fastq" 2> "$FILTERED_DIR/reports/${SRR}.txt"
        rm "$FASTQ" # liberar espacio
        CLEAN_FASTQ="$FILTERED_DIR/${SRR}_clean.fastq"
    else
        bowtie --quiet -v 3 -x "$RRNA_IDX" -q -1 "$FASTQ1" -2 "$FASTQ2" --un "$FILTERED_DIR/${SRR}_clean.fastq" 2> "$FILTERED_DIR/reports/${SRR}.txt"
        rm "$FASTQ1" "$FASTQ2"
        CLEAN_FASTQ="$FILTERED_DIR/${SRR}_clean.fastq"
    fi

    # 4. Mapeo con hisat2
    # Construir opciones de strandedness si corresponde
    STRAND_OPT=""
    if [ "$STRAND" != "NA" ]; then
        STRAND_OPT="--rna-strandness $STRAND"
    fi
    
        if [ "$TYPE" == "Single" ]; then
        # Single-end
        hisat2 --quiet $STRAND_OPT -x "$GENOMIC_IDX" -U "$CLEAN_FASTQ" -S "$HISAT_OUT/${SRR}.sam" 2> "$HISAT_OUT/reports/${SRR}.txt"
        rm "$CLEAN_FASTQ"
    else
        # Paired-end: se asume que tienes archivos _1 y _2 para cada SRR Paired-end
        CLEAN_FASTQ1="${FILTERED_DIR}/${SRR}_clean_1.fastq"
        CLEAN_FASTQ2="${FILTERED_DIR}/${SRR}_clean_2.fastq"
        if [ -f "$CLEAN_FASTQ1" ] && [ -f "$CLEAN_FASTQ2" ]; then
            hisat2 --quiet $STRAND_OPT -x "$GENOMIC_IDX" -1 "$CLEAN_FASTQ1" -2 "$CLEAN_FASTQ2" -S "$HISAT_OUT/${SRR}.sam" 2> "$HISAT_OUT/reports/${SRR}.txt"
            rm "$CLEAN_FASTQ1" "$CLEAN_FASTQ2"
        else
            echo "No se encontraron archivos Paired-end filtrados para $SRR"
            exit 1
        fi
    fi

    # 5. Convertir y ordenar .sam a .bam
    samtools view -b "$HISAT_OUT/${SRR}.sam" | samtools sort -o "$HISAT_OUT/${SRR}.bam"
    rm "$HISAT_OUT/${SRR}.sam"

    # 6. Cuantificar con Stringtie
    stringtie "$HISAT_OUT/${SRR}.bam" --rf -e -G "$GFF" -o "$STRINGTIE_OUT/${SRR}.gtf" -A "$STRINGTIE_OUT/gene_abundance/${SRR}.tab"

    # 7. Extraer FPKM y TPM
    cut -f1,8 "$STRINGTIE_OUT/gene_abundance/${SRR}.tab" | sort > "$STRINGTIE_OUT/gene_abundance/${SRR}_FPKM.txt"
    cut -f1,9 "$STRINGTIE_OUT/gene_abundance/${SRR}.tab" | sort > "$STRINGTIE_OUT/gene_abundance/${SRR}_TPM.txt"

    # Limpieza adicional si es necesario
    # rm "$HISAT_OUT/${SRR}.bam"
done < "$SAMPLES_FILE"

echo "Pipeline finalizado. Los resultados están en $STRINGTIE_OUT/gene_abundance/"