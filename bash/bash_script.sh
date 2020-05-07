#!/bin/sh

# installing the required packages  

apt-get update

apt install -y wget gcc libssl-dev make unzip openjdk-8-jdk git-all bwa samtools awscli virtualenv mysql-client mysql-server python3-dev

# install python (if not shipped with os) 
# cd /usr/src
# sudo wget https://www.python.org/ftp/python/3.6.6/Python-3.6.6.tar.xz
#sudo tar -xvf Python-3.6.6.tar.xz
#cd Python-3.6.6
#sudo ./configure --prefix=/opt/python3
#sudo make altinstall
#sudo ln -s /opt/python3/bin/python3.6 /usr/bin/python3.6

cd /opt 
mkdir software 
cd software 

# Installing GATK toolkit (Picard tools are bundled in with this) 
wget https://github.com/broadinstitute/gatk/releases/download/4.1.6.0/gatk-4.1.6.0.zip
unzip gatk-4.1.6.0.zip
rm -r gatk-4.1.6.0.zip

# Pulling the code into the instance 
cd /home/ubuntu
git clone https://github.com/Becky-Sadler/aws-ec2-pipeline.git
cd aws-ec2-pipeline
rm -r bash
mkdir run-directory reference
mv FH.bed reference 

# Sorting out reference sequence
cd reference

wget https://1000genomes.s3.amazonaws.com/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
gunzip hs37d5.fa.gz 

# creating .dict
/opt/software/gatk-4.1.6.0/gatk CreateSequenceDictionary -R hs37d5.fa

# creating .fa.fai
samtools faidx hs37d5.fa

# index for use with BWA

bwa index hs37d5.fa
