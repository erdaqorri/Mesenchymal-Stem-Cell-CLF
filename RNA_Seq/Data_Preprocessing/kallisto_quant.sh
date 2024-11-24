# Pseudo-align the reads to the reference Canis lupus familiaris transcriptome

READS_DIR="/path/to/fastp_output/"
INDEX="/path/to/Canis_lupus_familiaris.ROS_Cfam_1.0.cdna.all_index"
THREADS=10

if [ ! -d "${READS_DIR}" ]; then
    echo "Error: ${READS_DIR} does not exist"
    exit 1
fi


for R1_FILE in "${READS_DIR}"/*_R1_fastp.fastq.gz; do
    if [ ! -e "${R1_FILE}" ]; then
        echo "Error: ${R1_FILE} does not exist"
        exit 1
    elif [ ! -r "${R1_FILE}" ]; then
        echo "Error: ${R1_FILE} is not readable"
        exit 1
    fi

    SAMPLE=$(basename "${R1_FILE}" _R1_fastp.fastq.gz)

    R2_FILE="${READS_DIR}/${SAMPLE}_R2_fastp.fastq.gz"
    if [ ! -e "${R2_FILE}" ]; then
        echo "Error: ${R2_FILE} does not exist"
        exit 1
    elif [ ! -r "${R2_FILE}" ]; then
        echo "Error: ${R2_FILE} is not readable"
        exit 1
    fi

    OUTPUT_DIR="/path/to/kallisto_output/${SAMPLE}"

    if [ ! -d "${OUTPUT_DIR}" ]; then
        mkdir -p "${OUTPUT_DIR}"
    fi

    LOG_FILE="${OUTPUT_DIR}/${SAMPLE}.log"

    kallisto quant --index="${INDEX}" --output-dir="${OUTPUT_DIR}" --threads="${THREADS}" "${R1_FILE}" "${R2_FILE}" &> "${LOG_FILE}"
done
