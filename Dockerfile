FROM ubuntu:20.04
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y openjdk-8-jdk libncurses5-dev libbz2-dev lzma liblzma-dev zlib1g-dev \
    cmake unzip wget nano python3-pip python3
RUN wget -O bwa-0.7.17.tar.bz2 https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2/download \
    && tar xfv bwa-0.7.17.tar.bz2 && cd bwa-0.7.17 && make && ln -s /bwa-0.7.17/bwa /usr/local/bin
RUN wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 && \
    tar xfv samtools-1.10.tar.bz2 && rm samtools-1.10.tar.bz2 && cd samtools-1.10 && ./configure && make && make install
RUN wget -O gatk.zip "https://github.com/broadinstitute/gatk/releases/download/4.2.2.0/gatk-4.2.2.0.zip" && \
    unzip gatk.zip && rm gatk.zip && mv /gatk-4.2.2.0 /gatk
RUN pip install numpy biopython pandas matplotlib
ENV MPLCONFIGDIR=/tmp
RUN ln /usr/bin/python3 /usr/bin/python



RUN wget -O /opt/snpEff_latest_core.zip http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip/download && \
    cd /opt && unzip /opt/snpEff_latest_core.zip && rm -r clinEff snpEff_latest_core.zip
RUN pip install tqdm
COPY genes.gbk /opt/snpEff/data/covid19/genes.gbk

RUN echo "covid19.genome : covid19" >> /opt/snpEff/snpEff.config
RUN echo "  covid19.chromosomes : MN996528" >> /opt/snpEff/snpEff.config
RUN cd /opt/snpEff && java -jar snpEff.jar build -genbank -v covid19

RUN pip install scipy

COPY  minority_analysis.py /usr/local/bin/