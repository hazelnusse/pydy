import pexpect

def run_command(command):
    child.sendline(command)
    child.expect('(?<!-> )(\(\d+\))')
    return child.before + child.match.groups()[0]


child = pexpect.spawn('al')
child.expect('\(1\)')
print child.before
print child.after
print run_command('variables x,y,z')
print run_command('test = (x+y+z)^2')
print run_command('test2 = expand(test, 1:3)')
