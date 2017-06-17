from os import path


rule qdnaseq:
    input:
        expand('bam/deduped/{sample}.bam', sample=get_samples())
    output:
        logratios='qdnaseq/logratios.txt',
        calls='qdnaseq/calls.txt'
    params:
        bin_size=config['qdnaseq']['bin_size'],
        read_length=config['qdnaseq']['read_length'],
        options=config['qdnaseq']['extra'],
        spleens=get_spleen_samples()
    log:
        'logs/qdnaseq.log'
    run:
        output_dir = path.dirname(output.logratios)

        shell(
            'Rscript --vanilla scripts/qdnaseq.R'
            ' --out_dir {output_dir}'
            ' --bin_size {params.bin_size}'
            ' --read_length {params.read_length}'
            ' --spleens {params.spleens}'
            ' {params.options} {input} &> {log}')
