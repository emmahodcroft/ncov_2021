localrules: colors, download_for_cluster


ruleorder: finalize_swiss > finalize
ruleorder: filter_cluster > subsample
ruleorder: copy_from_scicore_archive > filter
ruleorder: rename_subclades_birds > rename_subclades
ruleorder: copy_from_scicore > download_for_cluster
ruleorder: copy_from_scicore_archive > download_for_cluster
ruleorder: copy_from_scicore_archive > copy_from_scicore
#ruleorder: download_masked > mask
#ruleorder: download_masked > diagnostic

def _get_path_for_cluster_input(cluster_wildcard):
    input_file = "{}/clusters/cluster_{}.txt".format(config.get("profile-name", ""), cluster_wildcard)

    if input_file:
        return path_or_url(input_file, keep_local=True)


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

rule extract_cluster_2:
    input:
        cluster = lambda wildcards: _get_path_for_cluster_input(wildcards.build_name), #"cluster_profile/clusters/cluster_{build_name}.txt",
        metadata = _get_unified_metadata,
        alignment = _get_unified_alignment,
        index = rules.index_sequences.output.sequence_index
    output:
        cluster_sample = "results/{build_name}/sample-precluster22.fasta"
    log:
        "logs/subsample_{build_name}_cluster-extract.txt"
    conda: config["conda_environment"]
    shell:
        """
        augur filter \
            --sequences {input.alignment} \
            --sequence-index {input.index} \
            --metadata {input.metadata} \
            --include {input.cluster} \
            --output {output.cluster_sample} 2>&1 | tee {log}
        """

rule filter_cluster:
    input:
        sequences = rules.extract_cluster.output.cluster_sample,
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
            --group-by country year month \
            --subsample-max-sequences 5000 \
            --output-strains {output.cluster} \
            --output {output.sequences} 2>&1 | tee {log}
        """
        #cp {input.cluster:q} {output.cluster:q}

rule copy_from_scicore:
    message: "copying files from Cornelius' runs"
    output:
        sequences = "results/filtered_gisaid.fasta.xz",
        metadata = "data/metadata.tsv",
        mutations = "results/mutation_summary_gisaid.tsv"
    conda: config["conda_environment"]
    shell:
        """
        cp ../../roemer0001/ncov-simple/pre-processed/gisaid/mutation_summary.tsv {output.mutations:q}
<<<<<<< HEAD
        cp ../../roemer0001/ncov-simple/data/gisaid/metadata.tsv {output.metadata:q}
        cp ../../roemer0001/ncov-simple/pre-processed/gisaid/filtered.fasta.xz {output.sequences:q}
=======
        cp ../../roemer0001/ncov-simple/pre-processed/metadata.tsv {output.metadata:q}
        xz -cdq ../../roemer0001/ncov-simple/pre-processed/gisaid/filtered.fasta.xz > {output.sequences:q}
>>>>>>> master
        """
        #cp ../../roemer0001/ncov-simple/data/gisaid/metadata.tsv {output.metadata:q}
        #xz -kz {output.mutations:q}
        #xz -kz {output.metadata:q}

rule copy_from_scicore_archive:
    message: "copying files from Cornelius' runs - the archive"
    output:
        sequences = "results/precomputed-filtered_gisaid.fasta",
        metadata = "data/downloaded_gisaid.tsv",
        mutations = "results/mutation_summary_gisaid.tsv"
    conda: config["conda_environment"]
    shell:
        """
        cp ../../roemer0001/ncov-simple/archive/pre-processed/gisaid/mutation_summary.tsv {output.mutations:q}
        cp ../../roemer0001/ncov-simple/archive/pre-processed/metadata.tsv {output.metadata:q}
        xz -cdq ../../roemer0001/ncov-simple/archive/pre-processed/gisaid/filtered.fasta.xz > {output.sequences:q}
        """
        #cp ../../roemer0001/ncov-simple/data/gisaid/metadata.tsv {output.metadata:q}
        #xz -kz {output.mutations:q}
        #xz -kz {output.metadata:q}
        #xz -cdq ../../roemer0001/ncov-simple/pre-processed/gisaid/filtered.fasta.xz > {output.sequences:q}

rule copy_from_scicore_archive:
    message: "copying files from Cornelius' runs - the archive"
    output:
        sequences = "results/filtered_gisaid.fasta.xz",
        metadata = "data/metadata.tsv",
        mutations = "results/mutation_summary_gisaid.tsv"
    conda: config["conda_environment"]
    shell:
        """
        cp ../../roemer0001/ncov-simple/archive/pre-processed/gisaid/mutation_summary.tsv {output.mutations:q}
        cp ../../roemer0001/ncov-simple/archive/pre-processed/metadata.tsv {output.metadata:q}
        cp ../../roemer0001/ncov-simple/archive/pre-processed/gisaid/filtered.fasta.xz {output.sequences:q}
        """
        #cp ../../roemer0001/ncov-simple/data/gisaid/metadata.tsv {output.metadata:q}
        #xz -kz {output.mutations:q}
        #xz -kz {output.metadata:q}

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
        #aws s3 cp s3://nextstrain-ncov-private/flagged-sequences_gisaid.tsv.xz - | xz -cdq > {output.flagged:q}
        #aws s3 cp s3://nextstrain-ncov-private/sequence-diagnostics_gisaid.tsv.xz - | xz -cdq > {output.diagnostics:q}
        #aws s3 cp s3://nextstrain-ncov-private/to-exclude.txt.xz - | xz -cdq > {output.to_exclude:q}
        #aws s3 cp s3://nextstrain-ncov-private/masked.fasta.xz - | xz -cdq > "results/masked.fasta"
        #aws s3 cp s3://nextstrain-ncov-private/masked.fasta.xz - | xz -cdq > {output.sequences:q}

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

