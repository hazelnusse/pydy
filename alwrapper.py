import re

import pexpect

class Al(object):

    def __init__(self):
        self._child = pexpect.spawn('al')
        self._child.expect('\(1\)')

    def run_command(self, command):
        """
        Runs the command in autolev and returns a list of (prompt, output).

        Example:
        >>>

        """
        self._child.sendline(command)
        self._child.expect(['(?<!-> )(\(\d+\))'])
        al_output = self._child.before
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


a = Al()
print a.run_command("variables a,b,c")
print a.run_command('test = (a+b+c)^2')
print a.run_command('test2 = expand(test, 1:3)')
print a.run_command('AUTOZ ON')
print a.run_command('frames D, E')
print a.run_command('simprot(E, D, 1, a)')
