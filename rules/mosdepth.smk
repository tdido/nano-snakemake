configfile: "config.yaml"

rule mosdepth_get:
    input:
        bam = "{aligner}/alignment/{sample}.bam",
        bai = "{aligner}/alignment/{sample}.bam.bai"
    threads: get_resource("mosdepth_get", "threads")
    resources:
        mem=get_resource("mosdepth_get", "mem"),
        walltime=get_resource("mosdepth_get", "walltime")
    output:
        protected("{aligner}/mosdepth/{sample}.mosdepth.global.dist.txt"),
        protected("{aligner}/mosdepth/{sample}.regions.bed.gz"),
    params:
        windowsize = 500,
        prefix = "{sample}",
        aligner = "{aligner}"
    log:
        "logs/{aligner}/mosdepth/mosdepth_{sample}.log"
    conda: "../envs/mosdepth.yaml"
    shell:
        "mosdepth --threads {threads} \
                  -n \
                  --by {params.windowsize} \
                  {params.aligner}/mosdepth/{params.prefix} {input.bam} 2> {log}"


rule mosdepth_combine:
    input:
        expand("{{aligner}}/mosdepth/{sample}.regions.bed.gz", sample=config["samples"])
    output:
        "{aligner}/mosdepth/regions.combined.gz"
    threads: get_resource("mosdepth_combine", "threads")
    resources:
        mem=get_resource("mosdepth_combine", "mem"),
        walltime=get_resource("mosdepth_combine", "walltime")
    log:
        "logs/{aligner}/mosdepth/mosdepth_combine.log"
    shell:
        os.path.join(workflow.basedir, "scripts/combine_mosdepth.py") + \
            " {input} -o {output} 2> {log}"


rule mosdepth_global_plot:
    input:
        expand("{{aligner}}/mosdepth/{sample}.mosdepth.global.dist.txt", sample=config["samples"])
    output:
        "{aligner}/mosdepth_global_plot/global.html"
    threads: get_resource("mosdepth_global_plot", "threads")
    resources:
        mem=get_resource("mosdepth_global_plot", "mem"),
        walltime=get_resource("mosdepth_global_plot", "walltime")
    log:
        "logs/{aligner}/mosdepth/mosdepth_global_plot.log"
    shell:
        os.path.join(workflow.basedir, "scripts/mosdepth_plot-dist.py") + \
            " {input} -o {output} 2> {log}"
