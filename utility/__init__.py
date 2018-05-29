import yaml
from os.path import expanduser, join


def load_config():
    home = expanduser("~")
    config_path = join(home, '.util.yaml')
    with open(config_path, 'r') as stream:
        try:
            my_config = yaml.load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    return my_config