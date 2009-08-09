from alwrapper import Autolev

def test_multiline():
    a = Autolev()
    assert a.run_command("variables a") == []
    assert a.run_command('AUTOZ ON') == []
    assert a.run_command('frames D, E') == []
    assert a.run_command('simprot(E, D, 1, a)') == [
            ('-> (5)', 'Z1 = COS(a)'),
            ('-> (6)', 'Z2 = SIN(a)'),
            ('-> (7)', 'E_D = [1, 0, 0; 0, Z1, -Z2; 0, Z2, Z1]')
            ]

def test_expand():
    s = "(-c3*s1 - c1*s2*s3)/((c1*c2 + c2*s1*(c3*s1 + c1*s2*s3)/(c1*c3 - s1*s2*s3))*(c1*c3 - s1*s2*s3))"
    s2 = s.replace("**", "^")

    a = Autolev()
    assert a.run_command("variables s1, s2, s3, c1, c2, c3") == []
    assert a.run_command('test = %s' % s2) == [
        ('-> (3)', 'test = -(c3*s1+c1*s2*s3)/(c2*c3*(c1^2+s1^2))')
        ]
    assert a.run_command('test2 = expand(test)') == [
            ('-> (5)',
                'test2 = -s1/(c2*(c1^2+s1^2)) - c1*s2*s3/(c2*c3*(c1^2+s1^2))')
            ]
