from os import path

import pandas as pd


rule gistic_prepare_markers:
    input:
        "qdnaseq/segmented.txt"
    params:
        chromosomes=config["references"]["gistic"]["chromosomes"],
        chromosome_map=config["references"]["gistic"]["chromosome_map"]
    output:
        "gistic/inputs/markers.txt"
    conda:
        path.join(workflow.basedir, "envs", "genopandas.yaml")
    script:
        path.join(workflow.basedir, "scripts", "gistic", "prepare_markers.py")


rule gistic_prepare_segments:
    input:
        "qdnaseq/segmented.txt"
    params:
        chromosomes=config["references"]["gistic"]["chromosomes"],
        chromosome_map=config["references"]["gistic"]["chromosome_map"]
    output:
        "gistic/inputs/segments.txt"
    conda:
        path.join(workflow.basedir, "envs", "genopandas.yaml")
    script:
        path.join(workflow.basedir, "scripts", "gistic", "prepare_segments.py")


rule gistic_prepare_samples:
    input:
        "gistic/inputs/segments.txt"
    params:
        samples=get_samples_for_group
    output:
        "gistic/inputs/samples.{group}.txt"
    run:
        # Check for missing samples (sanity check).
        segments = pd.read_csv(input[0], sep='\t')
        missing = set(params.samples) - set(segments['Sample'])

        if missing:
            raise ValueError('Missing samples {!r}'.format(missing))

        # Write sample list.
        samples = pd.DataFrame({'array': params.samples})
        samples.to_csv(output[0], sep='\t', index=False)


def _mcr_env(mcr_root=None, mcr_ver='v83'):
    if not mcr_root:
        return ""

    # LD path.
    ld_paths = ["{mcr_root}/{mcr_ver}/runtime/glnxa64",
                "{mcr_root}/{mcr_ver}/bin/glnxa64",
                "{mcr_root}/{mcr_ver}/sys/os/glnxa64",
                "${{LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}}"]

    ld_path = "LD_LIBRARY_PATH=" + ":".join(ld_paths)
    ld_path = ld_path.format(mcr_root=mcr_root, mcr_ver=mcr_ver)

    # XAPPLRESDIR path.
    xappl_path = ("XAPPLRESDIR={mcr_root}/{mcr_ver}"
                  "/MATLAB_Component_Runtime/{mcr_ver}/X11/app-defaults")
    xappl_path = xappl_path.format(mcr_root=mcr_root, mcr_ver=mcr_ver)

    return ld_path + " " + xappl_path


# Determine gistic and mcr root paths. If both are not given,
# we assume that gistic is available in PATH and that the
# LD_LIBRARY_PATH path has already been configured for the MCR.
gistic_root = config["rules"]["gistic"].get("gistic_root", "")
mcr_root = config["rules"]["gistic"].get("mcr_root", None)

if mcr_root is None and gistic_root:
    mcr_root = path.join(gistic_root, "MATLAB_Compiler_Runtime")

rule gistic:
    input:
        segments="gistic/inputs/segments.txt",
        markers="gistic/inputs/markers.txt",
        reference=config["references"]["gistic"]["mat"],
        samples="gistic/inputs/samples.{group}.txt"
    params:
        gistic_root=gistic_root,
        mcr_env=_mcr_env(mcr_root),
        output_dir="gistic/results/{group}",
        extra=" ".join(config["rules"]["gistic"]["extra"])
    output:
        amp_qplot="gistic/results/{group}/amp_qplot.pdf"
    log:
        "logs/gistic/{group}.log"
    shell:
        "{params.mcr_env} {params.gistic_root}gp_gistic2_from_seg"
        " -b {params.output_dir}"
        " -seg {input.segments}"
        " -refgene {input.reference}"
        " -alf {input.samples}"
        " {params.extra} > {log}"
