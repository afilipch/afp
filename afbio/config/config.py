# /usr/bin/python
'''Library provides API for configuration files loading, setting and parsing'''
import sys;
import os;

import yaml;

_confdir = os.path.dirname(os.path.realpath(__file__))

class LocalConfigError(Exception):
    pass;

#CONFIGS is used for shortcuts while calling configuration file
CONFIGS = {'samstat': 'samstat.yml', 'bedstat': 'bedstat.yml', 'chiflex': 'chiflex.yml', 'doublechiflex': 'doublechiflex.yml', 'chipchap': 'chipchap.yml', 'lrg': 'lrg.yml'}
for k, v in CONFIGS.items():
    CONFIGS[k] = os.path.join(_confdir, v)


def load_config(configuration):
    p = CONFIGS.get(configuration)
    if(not p):
        p = configuration;

    with open(p, 'r') as f:
        d = yaml.load(f)
    if(isinstance(d, dict)):
        return d;
    else:
        raise LocalConfigError('Config file %s is malformatted. It has to be convertible into dictionary' % path)



def save_config(obj, path):
    with open('path', 'w') as f:
            f.write(yaml.dump(obj, default_flow_style=False))


if (__name__=='__main__'):
    print(load(sys.argv[1]))