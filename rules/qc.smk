from os import path


rule multiqc:
    input:
        directory='.',
        fastqc=expand('qc/fastqc/{sample_lane}.{pair}_fastqc.html',
                      sample_lane=get_samples_with_lane(), pair=['R1']),
        samtools_stats=expand('qc/samtools_stats/{sample}.txt',
                              sample=get_samples())
    output:
        'qc/multiqc_report.html'
    params:
        ''
    log:
        'logs/multiqc.log'
    wrapper:
        'file://' + path.join(workflow.basedir, 'wrappers/multiqc')


rule fastqc:
    input:
        'fastq/trimmed/{sample}.{lane}.{pair}.fastq.gz'
    output:
        html='qc/fastqc/{sample}.{lane}.{pair}_fastqc.html',
        zip='qc/fastqc/{sample}.{lane}.{pair}_fastqc.zip'
    params:
        config['fastqc']['extra']
    wrapper:
        'file://' + path.join(workflow.basedir, 'wrappers/fastqc')


rule samtools_stats:
    input:
        'bam/deduped/{sample}.bam'
    output:
        'qc/samtools_stats/{sample}.txt'
    wrapper:
        'file://' + path.join(workflow.basedir, 'wrappers/samtools/stats')
