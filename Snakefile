segments = ['ha', 'na', 'pb2', 'pb1', 'pa', 'np', 'mp', 'ns']
lineages = ['h3n2', 'h1n1pdm']
resolutions = ['1y', '2y']

rule all:
    input:
        histogram_genome = expand('results/hist_distance_{lineage}_genome_{resolution}.png', lineage = lineages, resolution = '1y'),
        histogram_segments = expand('results/hist_distance_{lineage}_segments_{resolution}.png', lineage = lineages, resolution = '1y'),
        heatmap = expand('results/heatmap_{lineage}_{resolution}.png', lineage = lineages, resolution = '1y')

rule files:
    params:
        reference = 'config/reference_{lineage}_{segment}.gb'

files = rules.files.params

"""
rule align:
    message:
        '''
        Aligning sequences to {input.reference}
        - filling gaps with N
        '''
    input:
        sequences =
        reference = files.reference
    output:
        alignment = 'results/aligned_{lineage}_{segment}_{resolution}.fasta'
    shell:
        '''
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --fill-gaps \
            --remove-reference \
            --nthreads 1
        '''
"""

rule compare:
    message:
        '''
        Calculating pairwise genetic distance for all samples for segments and full-genome.
        '''
    input:
        alignments = expand('data/aligned_{{lineage}}_{segment}_{{resolution}}.fasta', segment = segments)
    output:
        pairwise = 'results/pairwise_{lineage}_genome_{resolution}.pickle'
    shell:
        '''
        python3 scripts/compare.py \
            --alignments {input.alignments} \
            --output {output.pairwise}
        '''

rule plot_histogram:
    message:
        '''
        Plotting histogram of genetic distance for segments & genomes.
        '''
    input:
        pairwise = rules.compare.output.pairwise
    params:
        cutoff = 30
    output:
        genome_histogram = 'results/hist_distance_{lineage}_genome_{resolution}.png',
        segments_histogram = 'results/hist_distance_{lineage}_segments_{resolution}.png'
    shell:
        '''
        python3 scripts/plot_hist.py \
            --pairwise {input.pairwise} \
            --cutoff {params.cutoff} \
            --lineage {wildcards.lineage} \
            --output-genome {output.genome_histogram} \
            --output-segments {output.segments_histogram}
        '''

rule plot_heatmap:
    message:
        '''
        Plotting heatmap of genetic distance for each segment.
        '''
    input:
        pairwise = rules.compare.output.pairwise
    params:
        max_distance = 300
    output:
        heatmap = 'results/heatmap_{lineage}_{resolution}.png'
    shell:
        '''
        python3 scripts/plot_heatmap.py \
            --pairwise {input.pairwise} \
            --max-distance {params.max_distance} \
            --lineage {wildcards.lineage} \
            --output {output.heatmap}
        '''
