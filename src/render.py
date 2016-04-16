#!/usr/bin/python
# -*- coding: utf-8 -*-

import cairo
from PIL import Image
import numpy as np
import math
from colour import Color


PI = np.pi
TWOPI = PI*2.
SIZE = 1000  # size of png image
NUM = 400  # number of nodes

BACK = 1.  # background color
GRAINS = 30
ALPHA = .3  # 0.05  # opacity of drawn points

ONE = 1./SIZE  # 1./SIZE


class Render(object):

    '''
    Inspired by @inconvergent (https://github.com/inconvergent),
    Render is a class to render images, for now it can draw points and lines,
    use ctx for all cairo context operations.
    '''

    def __init__(self):

        self._init_cairo()
        # self._get_colors(COLOR_PATH)
        self.colors = [(0, 0, 0)]
        self.n_colors = 1

    def _init_cairo(self):
        '''
        Initialize context and image
        '''
        sur = cairo.ImageSurface(cairo.FORMAT_ARGB32, SIZE, SIZE)
        ctx = cairo.Context(sur)
        ctx.scale(SIZE, SIZE)
        ctx.set_source_rgb(BACK, BACK, BACK)
        ctx.rectangle(0, 0, 1, 1)
        ctx.fill()

        self.sur = sur
        self.ctx = ctx

    def _get_colors(self, f):
        '''
        Misterious function, ask @inconvergent :)
        '''
        scale = 1./255.
        im = Image.open(f)
        w, h = im.size
        rgbim = im.convert('RGB')
        res = []
        for i in xrange(0, w):
            for j in xrange(0, h):
                r, g, b = rgbim.getpixel((i, j))
                res.append((r*scale, g*scale, b*scale))

        np.shuffle(res)
        self.colors = res
        self.n_colors = len(res)

    def draw_points(self, points):
        '''
        Draw a list of points
        '''
        r, g, b = self.colors[(NUM) % self.n_colors]
        self.ctx.set_source_rgba(255, 0, 0, ALPHA)
        for node in points:
            self.ctx.arc(node[0], node[1], 2*ONE,  -math.pi*2, 0)
            # self.ctx.rectangle(node[0], node[1], ONE, ONE)
            self.ctx.fill()

    def draw_lines(self, lines, color=None, filename=None):
        '''
        Draw a list of lines, black to red if color is True (FIXME),
        save a png for each line in filename if it isn't None
        '''
        for i in xrange(len(lines)):
            self.draw_line(lines[i][0], lines[i][1], color)
            if filename:
                self.sur.write_to_png(filename+str(i) + '.png')

    def draw_line(self, p1, p2, color=None):
        '''
        Draw a line from p1 to p2, black to red if color is True (FIXME)
        '''
        r, g, b = self.colors[(NUM) % self.n_colors]
        self.ctx.set_source_rgba(r, g, b, .5)
        steps = int(max([abs(p2[0] - p1[0]), abs(p2[1] - p1[1])]) / ONE)
        step1 = (max(p1[0], p2[0])-min(p1[0], p2[0]))/steps
        step1 = -step1 if p1[0] >= p2[0] else step1
        step2 = (max(p1[1], p2[1])-min(p1[1], p2[1]))/steps
        step2 = -step2 if p1[1] >= p2[1] else step2
        # colors = polylinear_gradient(['#0000ff','#ff0000'], steps+1)
        if color:
            #color = Color('black')
            colors = list(Color('black').range_to(Color('red'), steps+1))
            alpha = 1
        else:
            colors = list(Color('blue').range_to(Color('red'), steps+1))
            alpha = .2
        step = 0
        for x, y in zip(np.arange(p1[0], p2[0], step1), np.arange(p1[1], p2[1], step2)):
            self.ctx.set_source_rgba(*(colors[step].rgb + (alpha,)))
            # self.ctx.set_source_rgba(
            #    colors['r'][step], colors['g'][step], colors['b'][step], 1)
            step += 1
            self.ctx.rectangle(x, y, ONE, ONE)
            self.ctx.fill()
