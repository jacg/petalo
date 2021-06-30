from itertools import groupby, starmap
from operator import itemgetter
import argparse
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
import h5py
from sklearn.cluster import DBSCAN

parser = argparse.ArgumentParser()
parser.add_argument('infile', type=str, help='HDF5 file containing MC event data')
parser.add_argument('-n', '--min-cluster-size', type=int, default=10)
parser.add_argument('-d', '--max-distance-to_neighbour', type=float, default=100)
parser.add_argument('-q', '--charge-threshold', type=int, default=4, help='ignore sensors which detect fewer photons')
args = parser.parse_args()

MINIMUM_CLUSTER_SIZE          = args.min_cluster_size
MAXIMUM_DISTANCE_TO_NEIGHBOUR = args.max_distance_to_neighbour
CHARGE_THRESHOLD              = args.charge_threshold

def first_vertices_in_event(event_vertices):
    verts = tuple((evno, tid, x,y,z) for evno, tid, _, x,y,z, *crap, volume in event_vertices)
    def only(n): return lambda x: x[1]==n
    vsa = tuple(filter(only(1), verts))
    vsb = tuple(filter(only(2), verts))
    va = vsa[0][2:] if vsa else None
    vb = vsb[0][2:] if vsb else None
    if   va and vb: vs = np.array((va, vb))
    elif va       : vs = np.array((va,))
    elif vb       : vs = np.array((vb,))
    else          : vs = np.array
    return vs

def first_vertices(event_vertex_groups):
    for event_no, group in event_vertex_groups:
        yield event_no, first_vertices_in_event(group)

def synchronize_events(iterators, labels):
    while True:
        evnos_and_data = list(map(next, iterators))
        highest_event_number = max(map(itemgetter(0), evnos_and_data))
        for i, evn_and_data in enumerate(evnos_and_data):
            evn_and_data = list(evn_and_data) # enable mutation later in the loop
            while evn_and_data[0] < highest_event_number:
                if labels[i]: print(f'skipping event {evn_and_data[0]} in {labels[i]}')
                evn_and_data[0], evn_and_data[1] = next(iterators[i])
            evnos_and_data[i] = evn_and_data
        yield highest_event_number, tuple(map(itemgetter(1), evnos_and_data))


def cluster(hits, eps=MAXIMUM_DISTANCE_TO_NEIGHBOUR, min_samples=MINIMUM_CLUSTER_SIZE):
    label = DBSCAN(eps=eps, min_samples=min_samples).fit_predict(hits)
    return {n : hits[label==n] for n in set(label)}

def colour(label):
    if label ==  0: return 'yellow'
    if label ==  1: return 'pink'
    if label == -1: return 'gray'
    return 'gray'

def centroids_from_clusters(clusters):
    return np.array(tuple(h.mean(axis=0) for n,h in clusters.items() if n!=-1))


nema3_sources = np.array([[0,  10, 0], [0,  10, 56.25],
                          [0, 100, 0], [0, 100, 56.25],
                          [0, 200, 0], [0, 200, 56.25]])

def distance_between_point_and_line(point, line_end_1, line_end_2):
    line          = line_end_1 - line_end_2
    line_to_point =      point - line_end_2
    return norm(np.cross(line, line_to_point)) / norm(line)


def plotem(event_number, clusters, vers, active_source, miss):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    for n in sorted(clusters):
        h = clusters[n]
        ax.scatter(*h.T, c=colour(n), label=f'Cluster {n}' if n != -1 else 'noise')
    ax.scatter(*vers.T, marker='^', c='red', label='First vertex in LXe')
    centroids = centroids_from_clusters(clusters)
    #print(centroids)
    #print(vers)
    ax.scatter(*centroids.T, c='green', label='DBSCAN centroid')
    ax.set_title(f'event {event_number}: reco missed source by {miss:3.1f} mm')
    if len(vers)      == 2: ax.plot(*vers     .T, c='red')
    if len(centroids) == 2: ax.plot(*centroids.T, c='green')

    other_sources = np.array(tuple(filter(lambda p: not (p==active_source).all(), nema3_sources)))
    ax.scatter(*other_sources.T, c='palegreen', label='NEMA-3 sources')
    ax.scatter(*active_source  , c='red'      , label='LOR source')

    ax.legend()
    ax.set_xlim((-450,450)); ax.set_xlabel('x')
    ax.set_ylim((-450,450)); ax.set_ylabel('y')
    ax.set_zlim((-100,100)); ax.set_zlabel('z')
    plt.show()

def evhits(responses):
    for event_no, group in responses:
        yield event_no, np.array(tuple(sensor_xyz[sensor_id] for (event_no, sensor_id, charge) in group))



f = h5py.File(args.infile)
MC = f['MC']


sensor_xyz = {sensor: (x,y,z) for sensor, x,y,z in MC['sensor_xyz']}

primaries = iter(MC['primaries'])
vertices  = iter(MC['vertices'])
charges   = iter(MC['total_charge'])

significant_charges = filter(lambda x: x[2] >= CHARGE_THRESHOLD, charges)
event_sensor_responses = groupby(significant_charges, itemgetter(0))
event_hit_positions = evhits(event_sensor_responses)

event_vertex_groups = groupby(vertices, itemgetter(0))
event_first_vertices = first_vertices(event_vertex_groups)
event_primaries = ((evno, (x,y,z)) for evno, x,y,z, *p in primaries)

reco_misses = { tuple(pos): [] for pos in nema3_sources}

for (n,(v,h,p)) in synchronize_events((event_first_vertices, event_hit_positions, event_primaries),
                                       ("TRUE", "RECO", None)):
    c = cluster(h)
    centroids = centroids_from_clusters(c)
    if len(centroids) == 2 and len (v) == 2:
        active_source = np.array(p)
        reco_missed_by = distance_between_point_and_line(active_source, *centroids)
        true_missed_by = distance_between_point_and_line(active_source, *v)
        print(f'{true_missed_by:5.1f} {reco_missed_by:5.1f} {tuple(map(len,reco_misses.values()))}')
        reco_misses[p].append(reco_missed_by)
        plotem(n,c,v, active_source, reco_missed_by)
