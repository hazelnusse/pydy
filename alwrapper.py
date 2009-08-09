import re

import pexpect

def run_command(command):
    child.sendline(command)
    child.expect(['(?<!-> )(\(\d+\))'])
    al_output = child.before
    s = re.split("(-> \(\d+\))", al_output)
    input = s[0]
    assert input.strip() == command.strip()
    outputs = s[1:]
    pairs = []
    for i in range(len(outputs)/2):
        prompt = outputs[2*i]
        output = outputs[2*i+1]
        output = output.strip()
        pairs.append((prompt, output))
    return pairs
    #match = re.search("(-> \(\d+\))", al_output)
    #if match:
    #    input = al_output[:match.start()]
    #    prompt1 = match.groups()[0]
    #    output1 = al_output[match.end():]
    #else:
    #    input = al_output
    #    prompt1 = ""
    #    output1 = ""
    #return input, [prompt1, output1], child.match.groups()[0]


child = pexpect.spawn('al')
child.expect('\(1\)')
print child.before
print child.after
print run_command("variables a,b,c")
print run_command('test = (a+b+c)^2')
print run_command('test2 = expand(test, 1:3)')
print run_command('AUTOZ ON')
print run_command('frames D, E')
print run_command('simprot(E, D, 1, a)')
