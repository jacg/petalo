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


def test_synchronize_id():
    a = (1,2,3,4,  6,7,  9)
    b = (1,2,  4,5,6,  8,9)
    x = (1,2,  4,  6,    9)

    def ID(x): return x
    synced = tuple(synchronize((a,b), (ID, ID)))
    expected = tuple(((n, (n,n)) for n in x))
    print(synced)
    print(expected)
    assert synced == expected

def test_synchronize_attr():
    from operator import itemgetter
    a = ((1,10), (2,20),       (4,40),        (6,60))
    b = ((9,1 ),        (7,3), (6,4 ), (5,5), (4,6 ))
    expected = ((1, ((1,10), (9,1))),
                (4, ((4,40), (6,4))),
                (6, ((6,60), (4,6))))
    synced = tuple(synchronize((a,b), (itemgetter(0), itemgetter(1))))
    print(synced)
    print(expected)
    assert synced == expected
