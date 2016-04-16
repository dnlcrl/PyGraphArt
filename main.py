#!/usr/bin/python
# -*- coding: utf-8 -*-

from render import Render
from graph import Graph
import json
import numpy as np

FILENAME = 'output.png'


def load_italy_coords():
    '''load centered, normalized coords'''
    with open('data/comuni1.json') as comuni:
        data = json.load(comuni)

    coords = np.array([[x['lat'], x['lng']] for x in data], dtype=np.float32)
    return coords


def load_edges():
    with open('data/edges.json') as ef:
        edges = json.load(ef)
    return edges


def main():
    '''
    Create render and a graph, get some edge calling some graph
    function and draw them, finally save the results in FILENAME

    Bergamo = [.337125 ,.245148]
    Roma = [.4936765,.4637286]
    Napoli = [.5936468,.5253573]
    '''
    render = Render()
    #coords = load_italy_coords()
    g = Graph(
        points=None, oriented=False, rand=True, n=300, max_neighbours=9)

    edges = g.tsp(0)

    render.draw_points(g.nodes)
    render.draw_lines(g.coords_edges(g.edges))
    render.sur.write_to_png(FILENAME)
    import pdb
    pdb.set_trace()
    render.draw_lines(g.coords_edges(edges), color=True, filename='tsp/tsp')

    #render.draw_lines(g.coords_edges(edges), color=True, filename='kruskal/kru')
    render.sur.write_to_png(FILENAME)

if __name__ == '__main__':
    main()


def prepare_data():
    '''remove center and normalize values'''
    with open('data/comuni1.json') as comuni:
        data = json.load(comuni)

    coords = np.array([[x['lat'], x['lng']] for x in data], dtype=np.float32)
    angle = 90
    theta = (angle/180.) * np.pi
    rot_matrix = np.array(
        [[np.cos(theta), -np.sin(theta)], [np.sin(theta),  np.cos(theta)]])
    coords = np.dot(coords, rot_matrix)

    coords = coords - \
        np.array([np.max(coords, axis=0), np.min(coords, axis=0)]).mean(axis=0)
    coords /= 1.5*(np.max(coords, axis=0) - np.min(coords, axis=0))
    coords += [0.5, 0.5]

    for d, c in zip(data, coords):
        d['lat'] = c[0]
        d['lng'] = c[1]
    with open('data/comuni1.json', 'w') as comuni:
        json.dump(data, comuni)
    return coords
