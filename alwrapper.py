import os
import re

import pexpect

class AutolevError(Exception):
    pass

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
        if not input.strip() == command.strip():
            # autolev returned something more than just echoing the "command".
            # try if the echo was successful:
            if input.split("\n")[0].strip() == command.strip():
                # if yes, only raise an exception if something suspicious is
                # returned:
                if input.find("Error") != -1:
                    raise AutolevError("Unexpected result: %s" % input)
            else:
                # if not, always raise an exception
                raise AutolevError("Unexpected result: %s" % input)
        outputs = s[1:]
        pairs = []
        for i in range(len(outputs)/2):
            prompt = outputs[2*i]
            output = outputs[2*i+1]
            output = output.strip()
            pairs.append((prompt, output))
        return pairs
