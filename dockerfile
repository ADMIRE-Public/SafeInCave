# docker build -t safeincave .
# docker run --rm -it -v "C:\Users\herminiohon\Downloads\assignmentsafeincave-main:/app" safeincave bash

# Start from Ubuntu 22.04
FROM ubuntu:22.04

# Set noninteractive environment (avoid prompts in container)
ENV DEBIAN_FRONTEND=noninteractive

# Install Python3 and pip3
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    software-properties-common

# Set matplotlib backend to TkAgg
ENV MPLBACKEND=Agg

# Create a working directory
WORKDIR /app

# Copy your project into the container
COPY . /app

# Install Python dependencies using pip3
RUN pip3 install --upgrade pip
RUN pip3 install matplotlib
RUN pip3 install numpy==1.23.5
RUN pip3 install --no-cache-dir meshio==3.3.1
RUN pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu
RUN pip3 install pandas==1.4.3 
RUN apt install -y gmsh

# Add FEniCS 2019 PPA (Fenics 2019.1.0 example)
RUN apt-get install -y software-properties-common
RUN add-apt-repository ppa:fenics-packages/fenics
RUN apt-get update
RUN apt-get install -y fenics



