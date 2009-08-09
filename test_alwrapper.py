from alwrapper import Al

def test_multiline():
    a = Al()
    assert a.run_command("variables a") == []
    assert a.run_command('AUTOZ ON') == []
    assert a.run_command('frames D, E') == []
    assert a.run_command('simprot(E, D, 1, a)') == [
            ('-> (5)', 'Z1 = COS(a)'),
            ('-> (6)', 'Z2 = SIN(a)'),
            ('-> (7)', 'E_D = [1, 0, 0; 0, Z1, -Z2; 0, Z2, Z1]')
            ]
