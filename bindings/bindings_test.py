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

def test_implemented_in_petalo_fom():
    assert fulano.test_message() == "This is implemented in petalo::fom"
