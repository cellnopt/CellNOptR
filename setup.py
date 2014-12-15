import sys
import os
import easydev
# mimic setuptools -> develop, build, install
# build -> R CMD build
# develop or install -> R CMD install


def setup():
    argv = sys.argv
    # default mode is INSTALL
    if len(argv) == 1:
        mode = 'INSTALL'
    elif len(argv) == 2:
        easydev.check_param_in_list(argv[1], ['install', 'build'])
        mode = argv[1]

    for package in ['CellNOptR', 'CNORdt', 'CNORfuzzy', 'CNORfeeder', 'CNORode']:
        if mode == 'install':
            cmd = 'INSTALL'
        elif mode == 'build':
            cmd = 'build'
        cmd = 'R CMD %s packages' % cmd + os.sep + package
        print cmd
        os.system(cmd)

        
if __name__ == "__main__":
    setup()
