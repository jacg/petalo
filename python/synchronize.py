#!/usr/bin/env python3

from itertools import starmap
from operator import attrgetter

class Synchronizer:

    def __init__(self, iterable, key):
        self.data     = None
        self.count    = None
        self.key      = key
        self.iterator = iter(iterable)

    def step(self):
        self.data  = next(self.iterator)
        self.count = self.key(self.data)
        return self.count

def min_count(sources): return min(sources, key=attrgetter('count')).count
def max_count(sources): return max(sources, key=attrgetter('count')).count

def synchronize(iterables, keys):
    sources = tuple(starmap(Synchronizer, zip(iterables, keys)))
    while True:
        for source in sources:
            try                 : source.step()
            except StopIteration: return
        while min_count(sources) < (highest_count := max_count(sources)):
            for source in sources:
                if source.count < highest_count:
                     source.step()
        yield highest_count, tuple(map(attrgetter('data'), sources))
