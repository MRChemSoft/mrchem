from sys import float_info


def tuple_matches(t,
                  tolerance=1.0e-8,
                  error_definition='relative',
                  ignore_sign=False,
                  skip_below=float_info.min,
                  skip_above=float_info.max):
    """
    Checks if tuple matches based on tolerance settings.

    Returns:
        (tuple matches, error message) - error message is None if there is no error
    """

    assert error_definition in ['relative', 'absolute']

    x, x_ref = t

    # if ignore sign take absolute values
    if ignore_sign:
        x = abs(x)
        x_ref = abs(x_ref)

    if isinstance(x, int) and isinstance(x_ref, int):
        if x == x_ref:
            return (True, None)
        else:
            return (False, "expected: {0}".format(x_ref))

    if abs(x_ref) < skip_below:
        return (True, None)

    if abs(x_ref) > skip_above:
        return (True, None)

    error = x - x_ref
    if error_definition == 'relative':
        error /= x_ref

    if abs(error) <= tolerance:
        return (True, None)
    else:
        error_message = 'expected: {0} ({1} diff: {2:6.2e}'.format(x_ref, error_definition[:3], abs(error))
        if ignore_sign:
            error_message += ' ignoring signs'
        error_message += ')'
        return (False, error_message)


def test_tuple_matches():
    assert tuple_matches((13, 13)) == (True, None)
    assert tuple_matches((1.0 + 1.0e-9, 1.0)) == (True, None)
    assert tuple_matches((1.0 + 1.0e-9, 1.0), tolerance=1.0e-10) == (False, 'expected: 1.0 (rel diff: 1.00e-09)')
    assert tuple_matches((1.0 + 1.0e-7, 1.0)) == (False, 'expected: 1.0 (rel diff: 1.00e-07)')
    assert tuple_matches((0.01, 0.02), error_definition='absolute') == (False, 'expected: 0.02 (abs diff: 1.00e-02)')
    assert tuple_matches((0.01, 0.0002), error_definition='absolute', skip_below=0.001) == (True, None)
    assert tuple_matches((0.01, 2000.0), error_definition='absolute', skip_above=100.0) == (True, None)
    assert tuple_matches((10.0 + 1.0e-9, -10.0), error_definition='absolute') == (False, 'expected: -10.0 (abs diff: 2.00e+01)')
    assert tuple_matches((10.0 + 1.0e-9, -10.0), error_definition='absolute', ignore_sign=True) == (True, None)
    assert tuple_matches((13, 14)) == (False, 'expected: 14')
