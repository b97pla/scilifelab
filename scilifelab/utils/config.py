import ConfigParser
import os


def load_config(config_file=None):
    """Loads a configuration file.

    By default it assumes ~/.pm/pm.conf
    """
    try:
        if not config_file:
            config_file = os.path.join(os.environ.get('HOME'), '.pm', 'pm.conf')
        config = ConfigParser.SafeConfigParser()
        with open(config_file) as f:
            config.readfp(f)
        return config
    except IOError:
        raise IOError('There was a problem loading the configuration file. \
                Please make sure that ~/.pm/pm.conf exists and that you have \
                read permissions')
