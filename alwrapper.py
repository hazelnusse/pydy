import os
import re

import pexpect

class Autolev(object):

    def __init__(self):
        self._child = pexpect.spawn('al')
        self._child.expect('\(1\)')

    def __del__(self):
        # cleans the tmp output of autolev:
        os.system("rm -f alTmp.*")

    def run_command(self, command):
        """
        Runs the command in autolev and returns a list of (prompt, output).

        Example:
        >>> a = Al()
        >>> a.run_command("variables a")
        []
        >>> a.run_command('AUTOZ ON')
        []
        >>> a.run_command('frames D, E')
        []
        >>> a.run_command('simprot(E, D, 1, a)')
        [('-> (5)', 'Z1 = COS(a)'), ('-> (6)', 'Z2 = SIN(a)'), ('-> (7)', 'E_D
            = [1, 0, 0; 0, Z1, -Z2; 0, Z2, Z1]')]

        """
        self._child.sendline(command)
        self._child.expect(['(?<!-> )(\(\d+\))'])
        al_output = self._child.before
        s = re.split("(-> \(\d+\))", al_output)
        input = s[0]
        #print input
        #print command
        assert input.strip() == command.strip()
        outputs = s[1:]
        pairs = []
        for i in range(len(outputs)/2):
            prompt = outputs[2*i]
            output = outputs[2*i+1]
            output = output.strip()
            pairs.append((prompt, output))
        return pairs

#s = "(-c3*s1 - c1*s2*s3)/((c1*c2 + c2*s1*(c3*s1 + c1*s2*s3)/(c1*c3 - s1*s2*s3))*(c1*c3 - s1*s2*s3))"

#s2 = s.replace("**", "^")
##print s2
#a = Al()
#print a.run_command("variables a,b,c, s1, s2, s3, c1, c2, c3, u1, u2, u3, u4, u5, u6")
##print a.run_command('test = (a+b)^10')
#print a.run_command('test = %s' % s2)
#print a.run_command('test2 = expand(test)')
