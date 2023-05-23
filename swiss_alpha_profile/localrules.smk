localrules: colors, download_for_cluster


ruleorder: finalize_swiss > finalize
ruleorder: copy_base_files > extract_cluster
ruleorder: copy_base_files > filter_cluster
ruleorder: copy_base_files > subsample
ruleorder: copy_base_files > proximity_score
ruleorder: filter_cluster > subsample
ruleorder: rename_subclades_birds > rename_subclades

def _get_path_for_cluster_input(cluster_wildcard):
    input_file = "{}/clusters/cluster_{}.txt".format(config.get("profile-name", ""), cluster_wildcard)

    if input_file:
        return path_or_url(input_file, keep_local=True)

rule copy_base_files:
    message: "Copy the files that are unchanged between builds to save time"
    input:
        f1 = "results/starting_files/Alpha/sample-global.fasta",
        f2 = "results/starting_files/Alpha/sample-global.txt",
        f3 = "results/starting_files/Alpha/sample-precluster.fasta",
        f4 = "results/starting_files/Alpha/sample-cluster.fasta",
        f5 = "results/starting_files/Alpha/sample-cluster.txt",
        f6 = "results/starting_files/Alpha/proximity_cluster.tsv",
    output:
        f1 = "results/{build_name}/sample-global.fasta",
        f2 = "results/{build_name}/sample-global.txt",
        f3 = "results/{build_name}/sample-precluster.fasta",
        f4 = "results/{build_name}/sample-cluster.fasta",
        f5 = "results/{build_name}/sample-cluster.txt",
        f6 = "results/{build_name}/proximity_cluster.tsv",
    shell:
        """
        cp {input.f1} {output.f1}
        cp {input.f2} {output.f2}
        cp {input.f3} {output.f3}
        cp {input.f4} {output.f4}
        cp {input.f5} {output.f5}
        cp {input.f6} {output.f6}
        """

rule add_labels:
    message: "Remove extraneous colorings for main build and move frequencies"
    input:
        auspice_json = rules.incorporate_travel_history.output.auspice_json,
        tree = rules.refine.output.tree,
        clades = rules.clades.output.clade_data,
        mutations = rules.ancestral.output.node_data
    output:
        auspice_json = "results/{build_name}/ncov_with_accessions_and_travel_branches_and_labels.json",
    log:
        "logs/add_labels_{build_name}.txt"
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/add_labels.py \
            --input {input.auspice_json} \
            --tree {input.tree} \
            --mutations {input.mutations} \
            --clades {input.clades} \
            --output {output.auspice_json} 2>&1 | tee {log}
        """

rule finalize_swiss:
    message: "Remove extraneous colorings for main build and move frequencies"
    input:
        auspice_json = rules.add_labels.output.auspice_json,
        frequencies = rules.tip_frequencies.output.tip_frequencies_json,
        root_json = rules.export.output.root_sequence_json
    output:
        auspice_json = "auspice/ncov_{build_name}.json",
        root_json = "auspice/ncov_{build_name}_root-sequence.json",
        tip_frequency_json = "auspice/ncov_{build_name}_tip-frequencies.json"
    log:
        "logs/fix_colorings_{build_name}.txt"
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/fix-colorings.py \
            --input {input.auspice_json} \
            --output {output.auspice_json} 2>&1 | tee {log} &&
        cp {input.frequencies} {output.tip_frequency_json} &&
        cp {input.root_json} {output.root_json}
        """

rule extract_cluster:
    input:
        cluster = lambda wildcards: _get_path_for_cluster_input(wildcards.build_name), #"cluster_profile/clusters/cluster_{build_name}.txt",
        alignment = _get_unified_alignment
    output:
        cluster_sample = "results/{build_name}/sample-precluster.fasta"
    run:
        from Bio import SeqIO
        import lzma

        with open(input.cluster) as fh:
            cluster = set([x.strip() for x in fh.readlines()])

        seq_out = open(output.cluster_sample, 'w')
        with lzma.open(input.alignment, mode="rt") as fh:
            for s in SeqIO.parse(fh, 'fasta'):
                if s.id in cluster:
                    SeqIO.write(s, seq_out, 'fasta')

        seq_out.close()

rule filter_cluster:
    input:
        #sequences = rules.extract_cluster.output.cluster_sample,
        sequences = "results/{build_name}/sample-precluster.fasta",
        metadata = _get_unified_metadata,
        include = config["files"]["include"],
        cluster = lambda wildcards: _get_path_for_cluster_input(wildcards.build_name) #"{}/clusters/cluster_{build_name}.txt".format(config.get("profile-name",{})
    output:
        sequences = "results/{build_name}/sample-cluster.fasta",
        cluster = "results/{build_name}/sample-cluster.txt"
    log:
        "logs/subsample_{build_name}_cluster.txt"
    conda: config["conda_environment"]
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --include {input.include} \
            --output-strains {output.cluster} \
            --output {output.sequences} 2>&1 | tee {log}
        """
        #            --group-by year month \
        #    --subsample-max-sequences 9000 \
        # this is what is normally used if grouping by country, rather than just proximity
        # country grouping doesn't make sense here, because all are Swiss!
        #--group-by country year month \

rule download_for_cluster:
    message: "Downloading metadata and fasta files from S3"
    output:
        sequences = "results/precomputed-filtered_gisaid.fasta",
        metadata = "data/downloaded_gisaid.tsv",
        mutations = "results/mutation_summary_gisaid.tsv"
        #to_exclude = "results/to-exclude.txt"
        #diagnostics = "results/sequence-diagnostics_gisaid.tsv",
        #flagged = "results/flagged-sequences_gisaid.tsv",
    conda: config["conda_environment"]
    shell:
        """
        aws s3 cp s3://nextstrain-ncov-private/mutation-summary.tsv.xz - | xz -cdq > {output.mutations:q}
        aws s3 cp s3://nextstrain-ncov-private/metadata.tsv.gz - | gunzip -cq > {output.metadata:q}
        aws s3 cp s3://nextstrain-ncov-private/filtered.fasta.xz - | xz -cdq > {output.sequences:q}
        """

rule rename_subclades_birds:
    input:
        node_data = rules.subclades.output.clade_data
    output:
        clade_data = "results/{build_name}/subclades.json"
    run:
        import json
        with open(input.node_data, 'r', encoding='utf-8') as fh:
            d = json.load(fh)
            new_data = {}
            for k,v in d['nodes'].items():
                if "clade_membership" in v:
                    new_data[k] = {"Q677_membership": v["clade_membership"]}
        with open(output.clade_data, "w") as fh:
            json.dump({"nodes":new_data}, fh)