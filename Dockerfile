# Use Ubuntu 18.04 for better compatibility with older tools
FROM ubuntu:18.04

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Add CRAN repository for newer R version
RUN apt-get update && apt-get install -y \
    software-properties-common \
    dirmngr \
    wget \
    && wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | apt-key add - \
    && add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" \
    && apt-get update

# Install system dependencies (without r-base-dev initially)
RUN apt-get update && apt-get install -y \
    curl \
    git \
    unzip \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    python2.7 \
    python2.7-dev \
    python-pip \
    perl \
    openjdk-8-jdk \
    r-base \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Install r-base-dev separately with --fix-missing
RUN apt-get update && apt-get install -y --fix-missing r-base-dev || \
    (echo "r-base-dev installation failed, continuing without it" && \
     apt-get install -y r-recommended && \
     rm -rf /var/lib/apt/lists/*)

RUN apt-get update && apt-get install -y time

RUN pip install numpy

# Set JAVA_HOME
ENV JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64

# Install R packages with newer versions
RUN R -e "install.packages(c('ggplot2', 'ggthemes'), repos='https://cran.r-project.org/', dependencies=TRUE)"

# Create tools directory
RUN mkdir -p /opt/tools

# Install samtools 1.9
WORKDIR /opt/tools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar -xjf samtools-1.9.tar.bz2 && \
    cd samtools-1.9 && \
    ./configure && \
    make && \
    make install && \
    cd .. && rm -rf samtools-1.9*

# Install BWA 0.7.17
RUN wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 && \
    tar -xjf bwa-0.7.17.tar.bz2 && \
    cd bwa-0.7.17 && \
    make && \
    cp bwa /usr/local/bin/ && \
    cd .. && rm -rf bwa-0.7.17*

# Install freebayes - use a working approach
RUN apt-get update && apt-get install -y cmake libhts-dev && \
    wget https://github.com/freebayes/freebayes/releases/download/v1.3.6/freebayes-1.3.6-linux-amd64-static.gz && \
    gunzip freebayes-1.3.6-linux-amd64-static.gz && \
    chmod +x freebayes-1.3.6-linux-amd64-static && \
    mv freebayes-1.3.6-linux-amd64-static /usr/local/bin/freebayes

# Install GATK 4.1.7.0
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.1.7.0/gatk-4.1.7.0.zip && \
    unzip gatk-4.1.7.0.zip && \
    rm gatk-4.1.7.0.zip && \
    ln -s /opt/tools/gatk-4.1.7.0/gatk /usr/local/bin/gatk

# Install SnpEff 4.3t (specific version as in your working setup)
RUN cd /tmp && \
    wget https://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip && \
    unzip snpEff_latest_core.zip && \
    mv snpEff /opt/ && \
    rm snpEff_latest_core.zip

# Install Beagle 5.2
RUN wget https://faculty.washington.edu/browning/beagle/beagle.25Nov19.28d.jar -O /opt/tools/beagle.jar

# Clone NGS.PRSS1-2caller repository
WORKDIR /opt
RUN git clone https://github.com/Shuhua-Group/NGS.PRSS1-2caller.git

# Set up the NGS.PRSS1-2caller environment
WORKDIR /opt/NGS.PRSS1-2caller
RUN find . -name "*.sh" -exec chmod +x {} \; && \
    find . -name "*.pl" -exec chmod +x {} \; && \
    find . -name "*.py" -exec chmod +x {} \;

RUN cd /opt/NGS.PRSS1-2caller && \
    sed -i 's/chr_label = "no_chr"/chr_label = "add_chr"/g' remapping.py

# Fix the chromosome naming issue in the reference files - commented out since files already have correct naming
# RUN if [ -f "data/GRCh38.ALT_PRSS.fa" ]; then \
#     # Get the actual chromosome name from the FASTA file \
#     ACTUAL_CHR=$(grep "^>" data/GRCh38.ALT_PRSS.fa | head -1 | sed 's/>//g' | cut -d' ' -f1) && \
#     echo "Original chromosome name in FASTA: $ACTUAL_CHR" && \
#     # Keep a backup and rename chromosome to PRSS1_PRSS2 for consistency \
#     cp data/GRCh38.ALT_PRSS.fa data/GRCh38.ALT_PRSS.fa.backup && \
#     sed -i "s/^>.*/>PRSS1_PRSS2/" data/GRCh38.ALT_PRSS.fa && \
#     echo "Renamed chromosome to PRSS1_PRSS2 in FASTA file"; \
# fi

# Instead, just print the current file status
RUN echo "Checking NGS.PRSS1-2caller repository structure:" && \
    ls -la data/ && \
    echo "FASTA file header:" && \
    head -1 data/GRCh38.ALT_PRSS.fa

# Set up SnpEff database directory structure
RUN mkdir -p /opt/NGS.PRSS1-2caller/data/GRCh38_ALT_PRSS && \
    mkdir -p /opt/NGS.PRSS1-2caller/data/genomes

# Copy and fix the GFF file to ensure chromosome names match
RUN if [ -f "data/genes.gff.gz" ]; then \
    # Extract the original GFF file \
    gunzip -c data/genes.gff.gz > /opt/NGS.PRSS1-2caller/data/GRCh38_ALT_PRSS/genes.gff && \
    # Replace any chromosome names with PRSS1_PRSS2 \
    sed -i 's/^[^\t]*/PRSS1_PRSS2/' /opt/NGS.PRSS1-2caller/data/GRCh38_ALT_PRSS/genes.gff && \
    echo "Fixed chromosome names in GFF file to PRSS1_PRSS2"; \
else \
    echo "Warning: genes.gff.gz not found, will need to be provided at runtime"; \
fi

# Copy the FASTA sequence file for SnpEff - files already have correct chromosome names
RUN if [ -f "data/GRCh38.ALT_PRSS.fa" ]; then \
    cp data/GRCh38.ALT_PRSS.fa /opt/NGS.PRSS1-2caller/data/GRCh38_ALT_PRSS/sequences.fa && \
    cp data/GRCh38.ALT_PRSS.fa /opt/NGS.PRSS1-2caller/data/genomes/GRCh38_ALT_PRSS.fa && \
    echo "Copied FASTA files for SnpEff database"; \
fi

# Create SnpEff configuration
RUN echo "# SnpEff configuration for PRSS1-PRSS2 analysis" > /opt/snpEff/snpEff.config && \
    echo "data.dir = /opt/NGS.PRSS1-2caller/data" >> /opt/snpEff/snpEff.config && \
    echo "" >> /opt/snpEff/snpEff.config && \
    echo "# GRCh38_ALT_PRSS genome" >> /opt/snpEff/snpEff.config && \
    echo "GRCh38_ALT_PRSS.genome : PRSS1_PRSS2" >> /opt/snpEff/snpEff.config && \
    echo "GRCh38_ALT_PRSS.chromosomes : PRSS1_PRSS2" >> /opt/snpEff/snpEff.config && \
    echo "GRCh38_ALT_PRSS.PRSS1_PRSS2.codonTable : Standard" >> /opt/snpEff/snpEff.config

# Skip database build at build time - do it at runtime instead
RUN echo "SnpEff database will be built at runtime"

# Create parameter.txt with correct paths
RUN echo "# Software Path" > parameter.txt && \
    echo "PYTHON=/usr/bin/python2.7" >> parameter.txt && \
    echo "PERL=/usr/bin/perl" >> parameter.txt && \
    echo "JAVA=/usr/bin/java" >> parameter.txt && \
    echo "R=/usr/bin/Rscript" >> parameter.txt && \
    echo "SAMTOOLS=/usr/local/bin/samtools" >> parameter.txt && \
    echo "BWA=/usr/local/bin/bwa" >> parameter.txt && \
    echo "GATK=/opt/tools/gatk-4.1.7.0/gatk" >> parameter.txt && \
    echo "FREEBAYES=/usr/local/bin/freebayes" >> parameter.txt && \
    echo "BEAGLE=/opt/tools/beagle.jar" >> parameter.txt && \
    echo "SNPEFF=/opt/snpEff/snpEff.jar" >> parameter.txt && \
    echo "SNPSIFT=/opt/snpEff/SnpSift.jar" >> parameter.txt && \
    echo "# Resource Path" >> parameter.txt && \
    echo "ALT_REF=/opt/NGS.PRSS1-2caller/data/GRCh38.ALT_PRSS.fa" >> parameter.txt && \
    echo "REF=/opt/NGS.PRSS1-2caller/data/GRCh38.ALT_PRSS.fa" >> parameter.txt && \
    echo "hg38_REF=/opt/NGS.PRSS1-2caller/data/GRCh38.ALT_PRSS.fa" >> parameter.txt && \
    echo "hg19_REF=/opt/NGS.PRSS1-2caller/data/GRCh38.ALT_PRSS.fa" >> parameter.txt

# Test R installation and packages
RUN echo "Testing R installation..." && \
    R --version && \
    R -e "packageVersion('ggplot2')" && \
    R -e "packageVersion('ggthemes')" && \
    echo "R packages installed successfully"

# Set working directory
WORKDIR /data

# Add tools to PATH
ENV PATH="/opt/tools:/opt/snpEff:/usr/local/bin:${PATH}"

CMD ["/bin/bash"]
