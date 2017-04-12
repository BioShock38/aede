import numpy as np
import scipy.stats as stats
from functools import partial

class LocalAncestrySimulator(object):

    def __init__(self, gen_map, pipeline, name=None):
        self.name = name
        self.gen_map = gen_map
        self.pipeline = pipeline

    def run(self, *args, **kwargs):
        return self.pipeline(self.gen_map, *args, **kwargs)

    def __str__(self):
        return self.name

def morgan_to_index(pos_morgan, gen_map):
    """
    Convert an array of position in morgans into an array of index thanks to the
    genetic map
    """
    return[np.searchsorted(gen_map, p) for p in pos_morgan]

def jumps_builder(gen_map, lambd=6):
    max_cm = gen_map[-1]
    jumps = np.cumsum(stats.expon.rvs(scale=1/float(lambd), size=1000))
    return morgan_to_index(jumps[jumps < max_cm], gen_map)

jumps_6 = partial(jumps_builder, lambd=6)

def cluster_builder(jumps, gen_map, a_beta=0.8, b_beta=0.1):
    # alpha for one individu
    param_1 = -a_beta*(b_beta + a_beta**2 - a_beta) / (b_beta)
    param_2 = (b_beta + a_beta**2 - a_beta)*(a_beta - 1)/ (b_beta)
    alpha = stats.beta.rvs(param_1, param_2)
    clusters_by_region = stats.bernoulli.rvs(alpha, size=(len(jumps)+1)) 
    clusters_by_snp = np.repeat(clusters_by_region, np.ediff1d([0] + jumps + [len(gen_map)]))
    return clusters_by_snp

cluster_default = partial(cluster_builder, a_beta=0.8, b_beta=0.1)

def pipeline_hapmix(gen_map, pop_1, pop_2, jumps_builder, cluster_builder):
    h_adm = []
    jump_adm = []
    snp_cluster_adm = []
    for h_1, h_2 in zip(pop_1, pop_2):
        choices = [h_1, h_2]
        jumps = jumps_builder(gen_map)
        snp_cluster = cluster_builder(jumps, gen_map)
        h_adm.append(np.array([choices[cluster][j] for j, cluster in enumerate(snp_cluster)]))
        jumps_pos = np.zeros(len(h_1), dtype=np.bool_)
        #jumps_pos[jumps] = True
        jump_adm.append(jumps_pos)
        snp_cluster_adm.append(snp_cluster)

    return (np.array(h_adm), np.array(jump_adm), np.array(snp_cluster_adm))

default_pipeline_hapmix = partial(pipeline_hapmix, jumps_builder=jumps_6, cluster_builder=cluster_default)

def cluster_multiple_pops(gen_map, weights):
    alphas = np.random.dirichlet(weights)
    clusters = np.argmax(np.random.multinomial(1, alphas, 100000), axis=1)
    return clusters

def jumps_multiple_pops(gen_map, lambdas, clusters):
    max_cm = gen_map[-1]
    jumps = np.cumsum([stats.expon.rvs(scale=1/float(lambdas[c]), size=1) for c in clusters])
    return morgan_to_index(jumps[jumps < max_cm], gen_map)

def pipeline_multiple_pops(gen_map, l_pop, lambdas, weights):
    h_adm = []
    clusters_adm = []
    for choices in zip(*l_pop):
        clusters_all = cluster_multiple_pops(gen_map, weights)
        jumps = jumps_multiple_pops(gen_map, lambdas, clusters_all)
        nb_jumps = len(jumps)
        clusters_filtered = clusters_all[:nb_jumps+1]
        nb_to_repeat = np.ediff1d([0] + jumps + [len(gen_map)])
        print((clusters_filtered.shape, nb_to_repeat.shape))
        clusters = np.repeat(clusters_filtered, nb_to_repeat)
        print(clusters.shape)
        h_adm.append(np.array([choices[cluster][j] for j, cluster in enumerate(clusters)]))
        clusters_adm.append(clusters)

    return (np.array(h_adm), np.array(clusters_adm))
