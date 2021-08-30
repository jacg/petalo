from synchronize import synchronize

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
