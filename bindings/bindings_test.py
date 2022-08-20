from pytest import mark, param, raises

parametrize = mark.parametrize
xfail       = mark.xfail

import fulano

def test_fib():
    assert fulano.fib(10) == 89

def test_fab():
    assert fulano.fab(10) == 89

def test_lift():
    l = fulano.Lift(10)
    assert l.height == 10
    l.up(2)
    assert l.height == 12
    l.down(10)
    assert l.height ==  2

def test_fib_doc():
    assert fulano.fib.__doc__ == 'The naive, recursive fibonacci implementation'

def test_fab_doc():
    assert fulano.fab.__doc__ == 'The iterative fibonacci implementation'

def test_Lift_doc():
    assert fulano.Lift.__doc__ == "It's a Lift: it goes up and down"

from collections import namedtuple
C3d = namedtuple('C3d', 'x,y,z')
C2d = namedtuple('C2d', 'x,y')

@parametrize('arg, expected',
             ((1         , 'Int(1)'),
              ('abc'     , 'String("abc")'),
              ((1,2)     , 'IntTuple(1, 2)'),
              (('def', 3), 'StringTuple("def", 3)'),
              (C3d(1,2,3), 'Coordinates3d(1, 2, 3)'),
         param(C2d(1,2)  , 'Coordinates2d(1, 2)', marks=xfail(reason="IntTuple picks it up earlier")),
              (dict(a=1) , TypeError)
              ),
)
def test_rust_enum_parameter(arg, expected):
    if expected is not TypeError:
        result = fulano.rust_enum_parameter(arg)
        assert result == expected
    else:
        with raises(expected) as excinfo:
            fulano.rust_enum_parameter(arg)
        msg = str(excinfo.value)
        print(msg)
        assert msg == """argument 'e': failed to extract enum RustyEnum ('Int | String | IntTuple | StringIntTuple | Coordinates3d | Coordinates2d')
- variant Int (Int): 'dict' object cannot be interpreted as an integer
- variant String (String): 'dict' object cannot be converted to 'PyString'
- variant IntTuple (IntTuple): 'dict' object cannot be converted to 'PyTuple'
- variant StringIntTuple (StringIntTuple): 'dict' object cannot be converted to 'PyTuple'
- variant Coordinates3d (Coordinates3d): 'dict' object has no attribute 'x'
- variant Coordinates2d (Coordinates2d): 'dict' object has no attribute 'x'"""

cylinderX = namedtuple('cylinderX', 'y,z,r')
cylinderY = namedtuple('cylinderY', 'x,z,r')
cylinderZ = namedtuple('cylinderZ', 'x,y,r')
sphere    = namedtuple('sphere', 'x,y,z,r')

@parametrize('roi, expected',
             ((cylinderX(1,2,3), 'X 1 2 3'),
              (cylinderY(4,5,6), 'Y 4 5 6'),
              (cylinderZ(7,8,9), 'Z 7 8 9'),
              (sphere(2,4,6,8), 'S 2 4 6 8')
              ))
def test_rois(roi, expected):
    assert fulano.roi(roi) == expected


def test_crcs():
    from math import sin, cos, pi

    def polar(r, phi): return (r * cos(phi), r * sin(phi))

    step = pi / 6
    roi_from_centre = 50
    hot, cold, bg_radius, bg_activity = 4, 0, 4, 1

    rois = ((cylinderZ(*polar(roi_from_centre,  2*step),  4.0),  hot),
            (cylinderZ(*polar(roi_from_centre,  4*step),  6.5),  hot),
            (cylinderZ(*polar(roi_from_centre,  6*step),  8.5),  hot),
            (cylinderZ(*polar(roi_from_centre,  8*step), 11.0),  hot),
            (cylinderZ(*polar(roi_from_centre, 10*step), 14.0), cold),
            (cylinderZ(*polar(roi_from_centre, 12*step), 18.5), cold))

    bg_rois = (cylinderZ(*polar(roi_from_centre,  1*step), bg_radius),
               cylinderZ(*polar(roi_from_centre,  3*step), bg_radius),
               cylinderZ(*polar(roi_from_centre,  5*step), bg_radius),
               cylinderZ(*polar(roi_from_centre,  7*step), bg_radius),
               cylinderZ(*polar(roi_from_centre,  9*step), bg_radius),
               cylinderZ(*polar(roi_from_centre, 11*step), bg_radius))

    size   = (180,) * 3
    voxels =  (60,) * 3
    cfg = fulano.FomConfig(rois, bg_rois, bg_activity, voxels, size)
    crcs = cfg.crcs((1.0,) * (60*60*60))
    assert crcs == [0] * 6


def test_hmm():
    from fulano import fom_config
    rois = ((cylinderZ(1,2,3), 1),
            (cylinderY(3,2,4), 2))
    bg_rois = (cylinderX(9,8,7),
               sphere(1,2,3,4))
    # x = fom_config(rois, bg_rois, 1)
    # assert x == ''
