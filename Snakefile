segments = ['ha', 'na', 'pb2', 'pb1', 'pa', 'np', 'mp', 'ns']
lineages = ['h3n2', 'h1n1pdm']
resolutions = ['1y', '2y']

rule all:
    input:
        histogram_genome = expand('results/figs/hist_distance_{lineage}_genome_{resolution}.png', lineage = lineages, resolution = '1y'),
        histogram_segments = expand('results/figs/hist_distance_{lineage}_segments_{resolution}.png', lineage = lineages, resolution = '1y'),
        heatmap = expand('results/figs/heatmap_{lineage}_{resolution}.png', lineage = lineages, resolution = '1y'),
        hist_tmrca = expand('results/figs/hist_tmrca_{lineage}_{resolution}.png', lineage = lineages, resolution = '1y'),
        scatter_tmrca = expand('results/figs/scatter_tmrca_{lineage}_{resolution}.png', lineage = lineages, resolution = '1y'),
        pca_pairs = expand('results/figs/pca_pairs_{lineage}_{resolution}.png', lineage = lineages, resolution = '1y'),
        pca_strains = expand('results/figs/pca_strains_{lineage}_{resolution}.png', lineage = lineages, resolution = '1y'),

rule compare:
    message:
        '''
        For {wildcards.lineage}, calculating pairwise genetic distance for all samples for segments and full-genome.
        '''
    input:
        alignments = expand('data/aligned_{{lineage}}_{segment}_{{resolution}}.fasta', segment = segments)
    output:
        pairwise = 'results/pairwise_{lineage}_genome_{resolution}.hdf5'
    shell:
        '''
        python3 scripts/compare.py \
            --alignments {input.alignments} \
            --output {output.pairwise}
        '''

rule cluster:
    message:
        '''
        For {wildcards.lineage}, using SciPy's connected_components to cluster based on a genetic distance < than cutoff.
        '''
    input:
        pairwise = rules.compare.output.pairwise
    params:
        cutoff = 30  #When changing this number, also change the number in output as well.
    output:
        json = 'results/clusters_30_{lineage}_{resolution}.json'
    shell:
        '''
        python3 scripts/cluster.py  \
            --pairwise {input.pairwise} \
            --cutoff {params.cutoff} \
            --output {output.json}
        '''

rule plot_histogram:
    message:
        '''
        Plotting histogram of genetic distance for segments & genomes for {wildcards.lineage}
        '''
    input:
        pairwise = rules.compare.output.pairwise,
        clusters = rules.cluster.output.json
    params:
        cutoff = 30
    output:
        genome_histogram = 'results/figs/hist_distance_{lineage}_genome_{resolution}.png',
        segments_histogram = 'results/figs/hist_distance_{lineage}_segments_{resolution}.png'
    shell:
        '''
        python3 scripts/plot_hist.py \
            --pairwise {input.pairwise} \
            --cutoff {params.cutoff} \
            --lineage {wildcards.lineage} \
            --output-genome {output.genome_histogram} \
            --output-segments {output.segments_histogram}
        '''

def clock_rate(w):
    '''
    For a given lineage, this functions returns a list containig clock rate for each segment.
    This list is ordered according to segments.
    '''
    rate = {
        'h3n2': [0.00356, 0.00298, 0.00230, 0.00177, 0.00226, 0.00140, 0.00082, 0.00193],
         'h1n1pdm' : [0.00329, 0.00342, 0.00224, 0.00188, 0.00235, 0.00196, 0.00209, 0.00278]
    }
    return rate[(w.lineage)]

def segment_length(w):
    '''
    For a given lineage, this function returns a list containing length of each segment.
    The list is ordered according to segments.
    '''
    length = {
        'h3n2' : [1701, 1436, 2310, 2311, 2192, 1537, 999, 864],
        'h1n1pdm' : [1752, 1432, 2316, 2317, 2192, 1541, 1002, 865]
    }
    return length[(w.lineage)]

rule plot_heatmap:
    message:
        '''
        Plotting heatmap of genetic distance for each segment for {wildcards.lineage}
        '''
    input:
        pairwise = rules.compare.output.pairwise
    params:
        max_distance = 300,
        clock_rate = clock_rate,
        segment_length = segment_length
    output:
        heatmap = 'results/figs/heatmap_{lineage}_{resolution}.png'
    shell:
        '''
        python3 scripts/plot_heatmap.py \
            --pairwise {input.pairwise} \
            --max-distance {params.max_distance} \
            --clock-rate {params.clock_rate} \
            --segment-length {params.segment_length} \
            --lineage {wildcards.lineage} \
            --output {output.heatmap}
        '''

rule residuals:
    message:
        '''
        For {wildcards.lineage}, finding, plotting & saving cumulative residuals
        for pairwise TMRCA between all flu segments.
        '''
    input:
        pairwise = rules.compare.output.pairwise
    params:
        clock_rate = clock_rate,
        segment_length = segment_length
    output:
        scatterplot = 'results/figs/scatter_tmrca_{lineage}_{resolution}.png',
        histogram = 'results/figs/hist_tmrca_{lineage}_{resolution}.png'
    shell:
        '''
        python3 scripts/find_residuals.py \
            --pairwise {input.pairwise} \
            --clock-rate {params.clock_rate} \
            --segment-length {params.segment_length} \
            --lineage {wildcards.lineage} \
            --output-scatter {output.scatterplot} \
            --output-hist {output.histogram}
        '''

rule pca:
    message:
        '''
        For {wildcards.lineage}, performing PCA on strains & strain pairs.
        '''
    input:
        pairwise = rules.compare.output.pairwise,
        clusters = rules.cluster.output.json
    params:
        cutoff = 30,
        components = 8
    output:
        pca_pairs = 'results/figs/pca_pairs_{lineage}_{resolution}.png',
        pca_strains = 'results/figs/pca_strains_{lineage}_{resolution}.png'
    shell:
        '''
        python3 scripts/perform_pca.py  \
            --pairwise {input.pairwise} \
            --cutoff {params.cutoff} \
            --components {params.components} \
            --lineage {wildcards.lineage} \
            --output-pairs {output.pca_pairs} \
            --output-strains {output.pca_strains}
        '''
