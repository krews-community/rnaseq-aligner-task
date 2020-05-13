FROM openjdk:8-jre-slim as base
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    unzip \
    python-pip \
    python3 \
    python3-pip \
    liblzma-dev \
    libbz2-dev \
    zlib1g-dev \
    python \
     git

RUN mkdir /software
WORKDIR /software
ENV PATH="/software:${PATH}"


RUN wget https://github.com/alexdobin/STAR/archive/2.5.1b.tar.gz && tar -xzf 2.5.1b.tar.gz && \
    cd STAR-2.5.1b && make STAR && rm ../2.5.1b.tar.gz


ENV PATH="/software/STAR-2.5.1b/bin/Linux_x86_64:${PATH}"

RUN wget https://github.com/deweylab/RSEM/archive/v1.2.31.zip && \
    unzip v1.2.31.zip && rm v1.2.31.zip && \
    cd RSEM-1.2.31 && make

ENV PATH="/software/RSEM-1.2.31:${PATH}"

RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar xvjf samtools-1.9.tar.bz2 && cd samtools-1.9 && \
    ./configure --without-curses --disable-lzma --disable-bz2 && \
    make && make install && cd .. && \
    rm -r samtools-1.9 && rm samtools-1.9.tar.bz2


FROM openjdk:8-jdk-alpine as build
COPY . /src
WORKDIR /src

RUN ./gradlew clean shadowJar

FROM base
RUN mkdir /app
COPY --from=build /src/build/*.jar /app/star.jar
