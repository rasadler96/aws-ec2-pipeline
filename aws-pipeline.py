from ruffus import *
import os
import ruffus.cmdline as cmdline
import yaml
import sys
import re
import glob
import decimal
import configparser
import time

config = yaml.safe_load(open("aws-config.yml"))


#REFERENCE
reference_list = (config['reference'])
human_decoy = reference_list['human_decoy']

bed_list = (config['bed'])
fh_bed = bed_list['FH']

software_list = (config['software'])
gatk_picard = software_list['gatk_picard']
#picard_SortSam = software_list['picard_SortSam']
#mark_duplicates = software_list['mark_duplicates']
#gatk_jar = software_list['gatk']

parser = cmdline.get_argparse(description='Small pipeline for aws')

parser.add_argument("--input")
options = parser.parse_args()

run_directory = options.input

#check valid input dir provided
if os.path.isdir(run_directory):
    os.chdir(run_directory)
    log_dir = os.path.join(run_directory, "logs")
    if not os.path.exists(log_dir):
        os.mkdir(log_dir, 0o755)

    extension = ("*.fastq.gz")
    input_files = []
    input_files.extend(glob.glob(run_directory+"/"+extension))
    #print(input_files)

    #align with bwa mem
    @transform(input_files, formatter("(OMGL_\d+_U\d{3}).*_L001_R1_001.fastq.gz"), "{1[0]}_bwa.sam", ["{basename[0]}","logs/{1[0]}_bwa.log"])
    def align_with_bwa(input_file, output_file, extras):
        '''Align with BWA MEM'''

        fastq_1 = input_file
        fastq_2 = input_file.replace('L001_R1_001', 'L001_R2_001')
        reference = "/home/ubuntu/aws-ec2-pipeline/reference/hs37d5.fa"

        basename_cols = extras[0].split('_')
        sample_id = basename_cols[0]+"_"+basename_cols[1]+"_"+basename_cols[2]

        lane_id = sample_id

        read_group = "@RG\\tID:%s\\tPL:ILLUMINA\\tSM:%s\\tLB:library\\tPU:platform_unit" % (lane_id, sample_id)

        sam_file = output_file
        log_file = extras[1]

        #the command line statement we want to execute
        os.system("bwa mem -M -t 16 -R '%s' %s %s %s > %s 2>%s" % (read_group, reference, fastq_1, fastq_2, sam_file, log_file))

    @transform(align_with_bwa, formatter("(OMGL_\d+_U(\d{3}))_bwa.sam"), ["{1[0]}_bwa_sorted.bam", "logs/{1[0]}_indexBam.log"], ["{1[0]}_bwa.bam", "logs/{1[0]}_sam2bam.log", "{1[0]}_bwa_sorted.bam", "logs/{1[0]}_picardSort.log"])
    def index_bam(input_file, output_files, extras):
        '''Convert SAM to BAM, sort by coordinate and index'''

        index_log = output_files[1]
        bam, sam2bam_log, sorted_bam, picardSort_log = extras
        #the command line statement we want to execute
        os.system("samtools view -bS %s > %s 2>%s" % (input_file, bam, sam2bam_log))
        os.system("%s --java-options '-Xmx2g -Djava.io.tmpdir=/tmp' SortSam -I %s -O %s -SO coordinate --TMP_DIR /tmp &>%s" % (gatk_picard, bam, sorted_bam, picardSort_log))
        os.system("samtools index %s &>%s" % (extras[2], index_log))

    @transform(index_bam,
    #@transform(analysis_bam_list,
                formatter("(OMGL_\d+_U(\d{3}))_bwa_sorted.bam"),
                ["{1[0]}_dedup.bam", "logs/{1[0]}_picardDuplicates.log", "logs/{1[0]}_picardDupMetrics.log"]
                )
    def dedup_bam(input_file, output_files):
        dedup_bam, picard_log, picard_metrics = output_files
        os.system('%s --java-options "-Xmx8g -Djava.io.tmpdir=/tmp" MarkDuplicates -I %s -O %s --METRICS_FILE %s --TMP_DIR /tmp &>%s' % (gatk_picard, input_file[0], dedup_bam, picard_metrics, picard_log))

    #index duplicated marked BAM
    @follows(dedup_bam)
    @transform(dedup_bam, suffix("_dedup.bam"), "_dedup.bam.bai")
    def index_dedup_bam(input_files, output_file):
        bam_to_index = input_files[0]
        os.system('samtools index %s' % (bam_to_index))

    @transform(dedup_bam, formatter("(OMGL_\d+_U(\d{3}))_dedup.bam"), "{1[0]}.vcf", "logs/{1[0]}_gatk_hc.log")
    def call_variants(inputs, output_file, log_file):
        input_bam = inputs[0]
        '''Call variants with GATK HaplotypeCaller'''
        bed = fh_bed
        os.system('%s HaplotypeCaller -A Coverage -R %s -I %s -L %s -O %s &>%s' % (gatk_picard, human_decoy, input_bam, bed, output_file, log_file))

else:
    sys.stderr.write("The directory does not exist! Please check and try again\n")
    sys.exit()

cmdline.run(options)
