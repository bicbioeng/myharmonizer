FROM continuumio/miniconda3

WORKDIR /app

COPY ./myharmonizer .

RUN conda env create -f myharmonizer/myharmonizer.yml

# Make RUN commands use the new environment:
RUN echo "conda activate myharmonizer" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]

WORKDIR /app/myharmonizer
