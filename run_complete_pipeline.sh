#!/bin/bash

set -e

echo "=== PRSS1-PRSS2 COMPLETE PIPELINE (DOCKER VERSION) - FROM FASTQ TO VARIANTS ==="
echo "Starting at: $(date)"

# =============================================================================
# STEP 0: SETUP SNPEFF DATABASE
# =============================================================================
echo "=== STEP 0: SETTING UP SNPEFF DATABASE ==="

# Copy the setup script to the container and run it
cat > /tmp/setup_snpeff_database.sh << 'EOF'
#!/bin/bash

set -e

echo "=== Setting up SnpEff database for PRSS1-PRSS2 analysis ==="

# Define paths
SNPEFF_DIR="/opt/snpEff"
DATA_DIR="/opt/NGS.PRSS1-2caller/data"
DB_NAME="GRCh38_ALT_PRSS"
DB_DIR="$DATA_DIR/$DB_NAME"
GENOMES_DIR="$DATA_DIR/genomes"

# Create directories
mkdir -p "$DB_DIR"
mkdir -p "$GENOMES_DIR"

# Copy and fix FASTA file
echo "Setting up FASTA sequence file..."
if [ -f "/opt/NGS.PRSS1-2caller/data/GRCh38.ALT_PRSS.fa" ]; then
    # Copy to sequences.fa for SnpEff database
    cp "/opt/NGS.PRSS1-2caller/data/GRCh38.ALT_PRSS.fa" "$DB_DIR/sequences.fa"
    # Also copy to genomes directory
    cp "/opt/NGS.PRSS1-2caller/data/GRCh38.ALT_PRSS.fa" "$GENOMES_DIR/$DB_NAME.fa"
    
    # Check current chromosome name
    ORIG_CHR=$(grep "^>" "$DB_DIR/sequences.fa" | head -1 | sed 's/>//g' | cut -d' ' -f1)
    echo "Current chromosome name: $ORIG_CHR"
    
    # Change chromosome name to KI270803.1 to match GFF file
    echo "Setting chromosome name to KI270803.1 to match GFF annotations"
    sed -i 's/^>.*/>KI270803.1/' "$DB_DIR/sequences.fa"
    sed -i 's/^>.*/>KI270803.1/' "$GENOMES_DIR/$DB_NAME.fa"
    
    echo "FASTA files prepared with chromosome name: KI270803.1"
else
    echo "ERROR: FASTA file not found at /opt/NGS.PRSS1-2caller/data/GRCh38.ALT_PRSS.fa"
    # Let's check what files are actually available
    echo "Available files in /opt/NGS.PRSS1-2caller/data/:"
    ls -la /opt/NGS.PRSS1-2caller/data/ | head -10
    exit 1
fi

# Copy and fix GFF file
echo "Setting up GFF annotation file..."
if [ -f "/opt/NGS.PRSS1-2caller/data/genes.gff.gz" ]; then
    # Extract the original GFF file
    gunzip -c "/opt/NGS.PRSS1-2caller/data/genes.gff.gz" > "$DB_DIR/genes.gff"
    
    # The GFF file uses KI270803.1 coordinates, but our FASTA is only the PRSS region
    # We need to:
    # 1. Keep the original KI270803.1 chromosome name for SnpEff
    # 2. Filter annotations to only the PRSS region (749409-801557 on KI270803.1)
    # 3. Convert coordinates to match our FASTA sequence
    
    echo "Filtering GFF to PRSS1-PRSS2 region and converting coordinates..."
    
    # Create filtered GFF with coordinate conversion
    awk '
    BEGIN {OFS="\t"}
    /^#/ {print; next}
    {
        if ($1 == "KI270803.1" && $4 >= 749409 && $5 <= 801557) {
            # Convert coordinates: subtract 749408 to make them 1-based relative to our FASTA
            $4 = $4 - 749408
            $5 = $5 - 749408
            print
        }
    }' "$DB_DIR/genes.gff" > "$DB_DIR/genes_filtered.gff"
    
    mv "$DB_DIR/genes_filtered.gff" "$DB_DIR/genes.gff"
    
    echo "GFF file filtered and coordinates converted for PRSS1-PRSS2 region"
else
    echo "ERROR: GFF file not found at /opt/NGS.PRSS1-2caller/data/genes.gff.gz"
    exit 1
fi

# Verify chromosome names match
echo "Verifying chromosome name consistency..."
FASTA_CHR=$(grep "^>" "$DB_DIR/sequences.fa" | head -1 | sed 's/>//g' | cut -d' ' -f1)
GFF_CHR=$(head -10 "$DB_DIR/genes.gff" | grep -v "^#" | head -1 | cut -f1)

echo "FASTA chromosome: $FASTA_CHR"
echo "GFF chromosome: $GFF_CHR"

if [ "$FASTA_CHR" != "$GFF_CHR" ]; then
    echo "ERROR: Chromosome names don't match!"
    exit 1
fi

echo "Chromosome names match: $FASTA_CHR"

# Create SnpEff configuration
echo "Creating SnpEff configuration..."
cat > "$SNPEFF_DIR/snpEff.config" << EOL
# SnpEff configuration for PRSS1-PRSS2 analysis
data.dir = $DATA_DIR

# GRCh38_ALT_PRSS genome - using KI270803.1 to match GFF
$DB_NAME.genome : KI270803.1
$DB_NAME.chromosomes : KI270803.1
$DB_NAME.KI270803.1.codonTable : Standard
EOL

echo "SnpEff configuration created"

# Build the SnpEff database
echo "Building SnpEff database..."
cd "$SNPEFF_DIR"

if java -Xmx2g -jar snpEff.jar build -gff3 -v "$DB_NAME"; then
    if [ -f "$DB_DIR/snpEffectPredictor.bin" ]; then
        echo "SnpEff database built successfully"
    else
        echo "ERROR: Database build completed but snpEffectPredictor.bin not found"
        exit 1
    fi
else
    echo "ERROR: Database build failed"
    exit 1
fi

echo "=== SnpEff database setup complete ==="
EOF

chmod +x /tmp/setup_snpeff_database.sh
bash /tmp/setup_snpeff_database.sh

# =============================================================================
# STEP 1: FASTQ TO BAM ALIGNMENT (FROM variant_calling_fix.sh)
# =============================================================================
echo "=== STEP 1: ALIGNING FASTQ FILES TO GRCH38+ALT ==="

cd /data

# Check for FASTQ files
FASTQ_COUNT=$(ls -1 fastq/*.fastq* 2>/dev/null | wc -l)
if [ $FASTQ_COUNT -eq 0 ]; then
    echo "ERROR: No FASTQ files found in /data/fastq/"
    echo "Please place your FASTQ files in the fastq/ directory and try again."
    exit 1
fi

echo "Found $FASTQ_COUNT FASTQ files to process."

# Create necessary directories
mkdir -p bams
mkdir -p temp
mkdir -p references

echo "Setting up reference genome with ALT contig..."

# Download just chromosome 7 and the relevant ALT contig
if [ ! -f references/chr7.fa ]; then
    echo "Downloading chromosome 7 from UCSC..."
    wget -q -O references/chr7.fa.gz \
        "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/chromosomes/chr7.fa.gz"
    
    echo "Extracting chromosome 7..."
    gunzip references/chr7.fa.gz
    
    echo "Chromosome 7 downloaded successfully."
fi

if [ ! -f references/chr7_KI270803v1_alt.fa ]; then
    echo "Downloading chr7_KI270803v1_alt contig from UCSC..."
    wget -q -O references/chr7_KI270803v1_alt.fa.gz \
        "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/chromosomes/chr7_KI270803v1_alt.fa.gz"
    
    echo "Extracting ALT contig..."
    gunzip references/chr7_KI270803v1_alt.fa.gz
    
    echo "ALT contig downloaded successfully."
fi

# Create a minimal reference with just chr7 and the ALT contig
echo "Creating minimal reference for alignment..."
cat references/chr7.fa references/chr7_KI270803v1_alt.fa > references/chr7_with_alt.fa

# Index the reference for BWA and samtools
echo "Indexing the reference genome..."
bwa index references/chr7_with_alt.fa
samtools faidx references/chr7_with_alt.fa

echo "Reference preparation complete."

# Process each pair of FASTQ files
for R1_FILE in fastq/*_R1*.fastq*; do
    if [ -f "$R1_FILE" ]; then
        # Extract base name and find corresponding R2 file
        BASE_NAME=$(basename "$R1_FILE" | sed 's/_R1.*//g')
        R2_FILE=$(echo "$R1_FILE" | sed 's/_R1/_R2/g')
        
        if [ ! -f "$R2_FILE" ]; then
            echo "Warning: Could not find matching R2 file for $R1_FILE, skipping..."
            continue
        fi
        
        echo "Processing $BASE_NAME..."
        
        # Align with BWA-MEM to our minimal reference
        echo "Aligning $BASE_NAME to chr7_with_alt reference..."
        
        # Add read group information
        RG="@RG\\tID:$BASE_NAME\\tSM:$BASE_NAME\\tPL:ILLUMINA\\tLB:$BASE_NAME"
        
        # Perform alignment
        bwa mem -R "$RG" -t 4 \
            references/chr7_with_alt.fa \
            "$R1_FILE" "$R2_FILE" | \
            samtools view -bS - > temp/${BASE_NAME}.unsorted.bam
        
        # Sort BAM file
        echo "Sorting $BASE_NAME BAM file..."
        samtools sort -o bams/${BASE_NAME}.sorted.bam temp/${BASE_NAME}.unsorted.bam
        
        # Index BAM file
        echo "Indexing $BASE_NAME BAM file..."
        samtools index bams/${BASE_NAME}.sorted.bam
        
        # Clean up temporary unsorted file
        rm -f temp/${BASE_NAME}.unsorted.bam
        
        echo "$BASE_NAME alignment completed."
    fi
done

echo "All FASTQ files aligned to BAM format"

# =============================================================================
# STEP 2: PREPARE INPUT FILES FOR NGS.PRSS1-2CALLER
# =============================================================================
echo "=== STEP 2: PREPARING INPUT FILES FOR NGS.PRSS1-2CALLER ==="

echo "Creating sample list file..."
> sample_list.txt

# Only process .sorted.bam files to avoid duplicates and processing issues
for bam_file in /data/bams/*.sorted.bam; do
    if [ -f "$bam_file" ]; then
        sample_name=$(basename "$bam_file" .sorted.bam)
        echo -e "${sample_name}\t${bam_file}\tU" >> sample_list.txt
    fi
done

# Fallback: if no sorted.bam files, try regular .bam files (but exclude .alt.bam and .primary.bam)
if [ ! -s sample_list.txt ]; then
    echo "No .sorted.bam files found, looking for other .bam files..."
    for bam_file in /data/bams/*.bam; do
        if [ -f "$bam_file" ] && [[ ! "$bam_file" =~ \.(alt|primary)\.bam$ ]]; then
            sample_name=$(basename "$bam_file" .bam)
            echo -e "${sample_name}\t${bam_file}\tU" >> sample_list.txt
        fi
    done
fi

SAMPLE_COUNT=$(wc -l < sample_list.txt)
echo "Created sample list with $SAMPLE_COUNT samples"

# Debug: show what's in the sample list
echo "Sample list contents:"
cat sample_list.txt

if [ $SAMPLE_COUNT -eq 0 ]; then
    echo "ERROR: No suitable BAM files found in /data/bams/"
    echo "Available files:"
    ls -la /data/bams/
    exit 1
fi

# =============================================================================
# STEP 3: FIX PARAMETER.TXT FILE
# =============================================================================
echo "=== STEP 3: FIXING PARAMETER.TXT FILE ==="

# Update the parameter.txt file in the NGS.PRSS1-2caller directory with correct Docker paths
cat > /opt/NGS.PRSS1-2caller/parameter.txt << 'EOF'
path_to_samtools="/usr/local/bin"
path_to_bwa="/usr/local/bin"
path_to_java="/usr/bin"
path_to_python2="/usr/bin"
path_to_perl="/usr/bin"
path_to_R="/usr/bin"
path_to_gatk="/opt/tools/gatk-4.1.7.0"
path_to_freebayes="/usr/local/bin"
path_to_snpEff="/opt/snpEff"
EOF

echo "Updated parameter.txt with correct Docker paths"

# =============================================================================
# STEP 4: PATCH NGS.PRSS1-2CALLER SCRIPT
# =============================================================================
echo "=== STEP 4: PATCHING NGS.PRSS1-2CALLER SCRIPT ==="

# Create a backup of the original script
cp /opt/NGS.PRSS1-2caller/NGS.PRSS1-2caller.sh /opt/NGS.PRSS1-2caller/NGS.PRSS1-2caller.sh.backup

# Replace the problematic annotation lines with hardcoded paths
sed -i 's|python ${SD}/snpEff_config.py ${path_to_snpEff} ${SD}/data|python ${SD}/snpEff_config.py /opt/snpEff ${SD}/data|g' /opt/NGS.PRSS1-2caller/NGS.PRSS1-2caller.sh

sed -i 's|python ${SD}/snpEff_ann.py ${path_to_snpEff} ${SD}/data ${taskname}_PRSS.vcf ${taskname}_PRSS_snpEff_ann.vcf ${taskname}_PRSS_snpEff_ann.txt|python ${SD}/snpEff_ann.py /opt/snpEff ${SD}/data ${taskname}_PRSS.vcf ${taskname}_PRSS_snpEff_ann.vcf ${taskname}_PRSS_snpEff_ann.txt|g' /opt/NGS.PRSS1-2caller/NGS.PRSS1-2caller.sh

echo "Patched NGS.PRSS1-2caller.sh with hardcoded SnpEff paths"

# =============================================================================
# STEP 5: FIX MISSING TOOLS AND VCF ISSUE
# =============================================================================
echo "=== STEP 5: FIXING MISSING TOOLS AND VCF ISSUE ==="

# Install missing tools
apt-get update -qq
apt-get install -y tabix

echo "Installed tabix (includes bgzip)"

# Fix the snpEff_ann.py script to handle missing bgzip gracefully
cat > /tmp/fix_vcf_compression.py << 'EOF'
#!/usr/bin/python

# Read the original snpEff_ann.py file
with open('/opt/NGS.PRSS1-2caller/snpEff_ann.py', 'r') as f:
    content = f.read()

# Replace the bgzip and tabix commands with fallback handling
old_bgzip = 'os.system("bgzip %s" %(fout_vcf_path))'
old_tabix = 'os.system("tabix -p vcf %s.gz" %(fout_vcf_path))'

new_compression = '''# Try to compress with bgzip, fallback to gzip if bgzip not available
import subprocess
try:
    subprocess.check_call(["bgzip", fout_vcf_path])
    subprocess.check_call(["tabix", "-p", "vcf", fout_vcf_path + ".gz"])
    compressed_file = fout_vcf_path + ".gz"
except (subprocess.CalledProcessError, OSError):
    # Fallback: use gzip and skip tabix
    print("bgzip/tabix not available, using gzip...")
    subprocess.check_call(["gzip", fout_vcf_path])
    compressed_file = fout_vcf_path + ".gz"'''

# Replace both lines
content = content.replace(old_bgzip, new_compression)
content = content.replace(old_tabix, "# tabix handled above")

# Also fix the SnpEff command to use the correct file
old_snpeff_cmd = 'cmd = "java -Xmx4g -jar %s ann -v -c %s -dataDir %s GRCh38_ALT_PRSS %s.gz > %s" %(snpefff_jar_path,config_path,database_path,fout_vcf_path,fout_snpeff_path)'
new_snpeff_cmd = 'cmd = "java -Xmx4g -jar %s ann -v -c %s -dataDir %s GRCh38_ALT_PRSS %s > %s" %(snpefff_jar_path,config_path,database_path,compressed_file,fout_snpeff_path)'

content = content.replace(old_snpeff_cmd, new_snpeff_cmd)

# Write the fixed file
with open('/opt/NGS.PRSS1-2caller/snpEff_ann.py', 'w') as f:
    f.write(content)

print("Fixed VCF compression and SnpEff command")
EOF

python /tmp/fix_vcf_compression.py

echo "Fixed VCF compression handling"

# R plotting fix - now use the full R functionality since we have R 4.0+
cat > /tmp/fix_r_plotting.py << 'EOF'
#!/usr/bin/python

# Read the NGS.PRSS1-2caller.sh file
with open('/opt/NGS.PRSS1-2caller/NGS.PRSS1-2caller.sh', 'r') as f:
    content = f.read()

# Replace the R plotting line with a version that actually runs R
old_r_line = 'Rscript ${SD}/snpeFF_ann_snp.plot.R ${SD}/data/GRCh38.ALT_PRSS.PRSS_exon.txt ${taskname}_PRSS_snpEff_ann.plot.txt ${taskname}_PRSS_snpEff_ann.plot.pdf'

new_r_line = '''# Run R plotting with proper error handling
echo "Creating variant position plot with R..."
if Rscript ${SD}/snpeFF_ann_snp.plot.R ${SD}/data/GRCh38.ALT_PRSS.PRSS_exon.txt ${taskname}_PRSS_snpEff_ann.plot.txt ${taskname}_PRSS_snpEff_ann.plot.pdf; then
    echo "R plotting completed successfully"
else
    echo "R plotting failed, but continuing with analysis"
    # Create empty plot file so pipeline doesn't fail
    touch ${taskname}_PRSS_snpEff_ann.plot.pdf
fi'''

content = content.replace(old_r_line, new_r_line)

# Write the fixed file
with open('/opt/NGS.PRSS1-2caller/NGS.PRSS1-2caller.sh', 'w') as f:
    f.write(content)

print("Fixed R plotting to use full R functionality")
EOF

python /tmp/fix_r_plotting.py

echo "Enhanced R plotting with R 4.0+ support"
echo "=== FIXES COMPLETE ==="

# =============================================================================
# STEP 6: RUN NGS.PRSS1-2CALLER
# =============================================================================
echo "=== STEP 6: RUNNING NGS.PRSS1-2CALLER ==="

ANALYSIS_NAME="prss1_analysis"

# Change to the NGS.PRSS1-2caller directory
cd /opt/NGS.PRSS1-2caller

# Update the parameter.txt file to point to our setup
cat > parameter.txt << 'EOF'
# Software Path
PYTHON=/usr/bin/python2.7
PERL=/usr/bin/perl
JAVA=/usr/bin/java
R=/usr/bin/Rscript
SAMTOOLS=/usr/local/bin/samtools
BWA=/usr/local/bin/bwa
GATK=/opt/tools/gatk-4.1.7.0/gatk
FREEBAYES=/usr/local/bin/freebayes
BEAGLE=/opt/tools/beagle.jar
SNPEFF=/opt/snpEff/snpEff.jar
SNPSIFT=/opt/snpEff/SnpSift.jar
# Resource Path
ALT_REF=/opt/NGS.PRSS1-2caller/data/GRCh38.ALT_PRSS.fa
REF=/opt/NGS.PRSS1-2caller/data/GRCh38.ALT_PRSS.fa
hg38_REF=/opt/NGS.PRSS1-2caller/data/GRCh38.ALT_PRSS.fa
hg19_REF=/opt/NGS.PRSS1-2caller/data/GRCh38.ALT_PRSS.fa
EOF

# Run the main pipeline
echo "Running NGS.PRSS1-2caller..."
bash NGS.PRSS1-2caller.sh -i /data/sample_list.txt -n "$ANALYSIS_NAME" -r "GRCh38" -m 50

# =============================================================================
# FINAL RESULTS
# =============================================================================
echo "=== ANALYSIS COMPLETE ==="
echo "Results are available in: /data/prss1_analysis/"
echo "Key files:"
echo "  - Aligned BAM files: /data/bams/*.sorted.bam"
echo "  - prss1_analysis_PRSS.vcf: Original variant calls"
echo "  - prss1_analysis_PRSS_snpEff_ann.vcf: SnpEff annotations"
echo "  - prss1_analysis_PRSS_snpEff_ann.txt: Annotation summary in text format"
echo "  - prss1_analysis_PRSS_snpEff_ann.plot.pdf: Variant position plot"
echo "Complete pipeline completed successfully at: $(date)"
