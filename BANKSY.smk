rule BANKSY_clustering:
    output:
        os.path.join("results","{ds}","BANKSY","{sp}","clustering.tsv"),
    params:
        script = os.path.join("workflow","scripts","BANKSY_run.py"),
        n = lambda wildcards: samples_df[(samples_df.dataset_id == wildcards.ds)
        & (samples_df.sample_id == wildcards.sp)].n_cluster.values[0],

    log:    
        os.path.join("results","{ds}","BANKSY","{sample}","BANKSY.log")
    conda:
        "../envs/BANKSY.yml"
    shell:
        ’python {params.script} {wildcards.ds} {wildcards.sp} {params.n}’
        ’&> {log}’