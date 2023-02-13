import pytest
import os
import shlex
import subprocess
import sys
import tempfile


MODULE = "import_serial"


def run(line, **kwargs):
    # try:
    #      if os.environ["TRAVIS"] == "true":
    #         # line += (" -t")
    #         # print("Light-testing for TRAVIS CI - not executing REFMAC5")
    #         pass
    # except KeyError:
    #     # print("Normal testing mode.")
    #     pass
    if sys.platform == 'win32':
        posix_bool = False
        shell_bool = True
    else:
        posix_bool = True
        shell_bool = False
    print('\n$ ccp4-python -m ' + MODULE, line)
    # print('\n$ python -m', MODULE, line)
    command = ["ccp4-python", '-m', MODULE] + shlex.split(line, posix=posix_bool)
    # return subprocess.run(command,
    process = subprocess.Popen(command,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               universal_newlines=True,
                               shell=shell_bool,
                               **kwargs)
    cp = AttrDict()
    cp.stdout, cp.stderr = process.communicate()
    cp.returncode = process.returncode
    return cp


def run_ok(*args, **kwargs):
    cp = run(*args, **kwargs)
    assert cp.returncode == 0
    assert not cp.stderr
    print(cp.stdout)
    return cp


#def config(name):
#    return os.path.dirname(os.path.realpath(__file__)) + "/fixtures/" + name
#    # return pathlib.Path(__file__).parent / 'fixtures' / name


@pytest.fixture
def tmp_environ():
    os.chdir(tempfile.gettempdir())
    return


# https://stackoverflow.com/
# questions/4984647/accessing-dict-keys-like-an-attribute
class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self
