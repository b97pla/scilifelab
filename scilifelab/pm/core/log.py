"""pmtools log core module"""

import os
import logbook
from logbook import Logger

from cement.core import log
from cement.utils.misc import is_true
from cement.ext.ext_logging import LoggingLogHandler
from cement.utils import fs

class PmLogHandler(log.CementLogHandler):  
    """
    PmLogHandler - override CementLogHandler to use logbook.

    This class uses the same configuration options as 
    :ref:`LoggingLogHandler <cement.ext.ext_logging>`
    """
    
    class Meta:
        interface = log.ILog
        """The interface that this class implements."""

        label = 'pmlog'
        """The string identifier of this handler."""

        namespace = "pm"
        """
        The logging namespace.  
        
        Note: Although Meta.namespace defaults to None, Cement will set this 
        to the application label (CementApp.Meta.label) if not set during
        setup.
        """

        file_format = "{record.time} ({record.level_name}) {record.channel} : {record.message}"
        """The logging format for the file logger."""

        console_format = "{record.time:%Y-%m-%d %H:%M} ({record.level_name}): {record.message}"
        """The logging format for the consoler logger."""

        debug_format = "{record.time} ({record.level_name}) {record.channel} : {record.message}"
        """The logging format for both file and console if ``debug==True``."""

        log_setup = None
        """Nested log setup placeholder"""

        level = 0
        """Global level for handlers"""

        clear_loggers = True
        """Whether of not to clear previous loggers first."""
        # These are the default config values, overridden by any '[log]' 
        # section in parsed config files.
        config_section = 'log'
        """
        The section of the application configuration that holds this handlers
        configuration.
        """
        
        config_defaults = dict(
            file=None,
            level='INFO',
            to_console=True,
            rotate=False,
            max_bytes=512000,
            max_files=4,
            )
        """
        The default configuration dictionary to populate the ``log`` section.
        """
            
    levels = ['INFO', 'WARN', 'ERROR', 'DEBUG', 'FATAL']


    def __init__(self, *args, **kw):
        super(PmLogHandler, self).__init__(*args, **kw)
        self.app = None
        
    def _setup(self, app_obj):
        super(PmLogHandler, self)._setup(app_obj)
        if self._meta.namespace is None:
            self._meta.namespace = self.app._meta.label

        self.backend = Logger(self._meta.namespace)

        # hack for application debugging
        if is_true(self.app._meta.debug):
            self.app.config.set('log', 'level', 'DEBUG')
            
        # Mainly for backwards compatibility since Logger level should
        # be NOTSET (level 0). Output level is controlled by handlers
        self.set_level(self.app.config.get('log', 'level'))
        
        # clear loggers?
        if is_true(self._meta.clear_loggers):
            self.clear_loggers()
            
        # console
        if is_true(self.app.config.get('log', 'to_console')):
            self._setup_console_log()
        
        # file
        if self.app.config.get('log', 'file'):
            self._setup_file_log()
        # nested setup
        self.backend.handlers.append(logbook.NullHandler(bubble=False))
        self.log_setup = logbook.NestedSetup(self.backend.handlers)
        with self._console_handler.applicationbound():
            self.debug("logging initialized for '%s' using PmLogHandler" % \
                           self._meta.namespace)


    def set_level(self, level):
        """
        Set the log level.  Must be one of the log levels configured in 
        self.levels which are ``['INFO', 'WARN', 'ERROR', 'DEBUG', 'FATAL']``.
        
        :param level: The log level to set.
        
        """
        level = level.upper()
        if level not in self.levels:
            level = 'INFO'
        level = logbook.lookup_level(level.upper())
        self.level = level
        
    def get_level(self):
        """Returns a string representation of the current log level."""
        return logbook.get_level_name(self.level)

    def _setup_console_log(self):
        """Add a console log handler."""
        if logbook.lookup_level(self.get_level()) == logbook.DEBUG:
            fmt_string = self._meta.debug_format
        else:
            fmt_string = self._meta.console_format
        console_handler = logbook.StderrHandler(
            format_string=fmt_string,
            level = logbook.lookup_level(self.get_level()),
            bubble = True)
        self._console_handler = console_handler
        self.backend.handlers.append(console_handler)

    def _setup_file_log(self):
        """Add a file log handler."""
        file_path = os.path.expandvars(fs.abspath(self.app.config.get('log', 'file')))
        log_dir = os.path.dirname(file_path)
        if not os.path.exists(log_dir):
            os.makedirs(log_dir)
        if logbook.lookup_level(self.get_level()) == logbook.DEBUG:
            fmt_string = self._meta.debug_format
        else:
            fmt_string = self._meta.file_format

        if self.app.config.get('log', 'rotate'):
            from logbook import RotatingFileHandler
            file_handler = RotatingFileHandler(
                file_path, 
                max_size=int(self.app.config.get('log', 'max_bytes')), 
                backup_count=int(self.app.config.get('log', 'max_files')),
                format_string=fmt_string,
                level = logbook.lookup_level(self.get_level()),
                bubble = True,
                )
        else:
            from logbook import FileHandler
            file_handler = FileHandler(file_path,
                                       format_string=fmt_string,
                                       level = logbook.lookup_level(self.get_level()),
                                       bubble = True,
                                       )
        
        self._file_handler = file_handler
        self.backend.handlers.append(file_handler)

    def _get_logging_kwargs(self, namespace, **kw):
        if namespace is None:
            namespace = self._meta.namespace
        if 'extra' in kw.keys() and 'namespace' in kw['extra'].keys():
            pass
        elif 'extra' in kw.keys() and 'namespace' not in kw['extra'].keys():
            kw['extra']['namespace'] = namespace
        else:
            kw['extra'] = dict(namespace=namespace)
        
        return kw

    def info(self, msg, namespace=None, **kw):
        """
        Log to the INFO facility.
        
        :param msg: The message the log.
        :param namespace: A log prefix, generally the module ``__name__`` that 
            the log is coming from.  Will default to self._meta.namespace if 
            None is passed.
        :keyword kw: Keyword arguments are passed on to the backend logging 
            system.
            
        """
        kwargs = self._get_logging_kwargs(namespace, **kw)
        self.backend.info(msg, **kwargs)

    def debug(self, msg, namespace=None, **kw):
        """
        Log to the DEBUG facility.
        
        :param msg: The message the log.
        :param namespace: A log prefix, generally the module ``__name__`` that 
            the log is coming from.  Will default to self._meta.namespace if 
            None is passed.  For debugging, it can be useful to set this to 
            ``__file__``, though ``__name__`` is much less verbose.
        :keyword kw: Keyword arguments are passed on to the backend logging 
            system.
                        
        """
        kwargs = self._get_logging_kwargs(namespace, **kw)
        self.backend.debug(msg, **kwargs)

    def warn(self, msg, namespace=None, **kw):
        """
        Log to the WARN facility.
        
        :param msg: The message the log.
        :param namespace: A log prefix, generally the module ``__name__`` that 
            the log is coming from.  Will default to self._meta.namespace if 
            None is passed.  For debugging, it can be useful to set this to 
            ``__file__``, though ``__name__`` is much less verbose.
        :keyword kw: Keyword arguments are passed on to the backend logging 
            system.
                        
        """
        kwargs = self._get_logging_kwargs(namespace, **kw)
        self.backend.warn(msg, **kwargs)

    def critical(self, msg, namespace=None, **kw):
        """
        Log to the CRITICAL facility.
        
        :param msg: The message the log.
        :param namespace: A log prefix, generally the module ``__name__`` that 
            the log is coming from.  Will default to self._meta.namespace if 
            None is passed.  For debugging, it can be useful to set this to 
            ``__file__``, though ``__name__`` is much less verbose.
        :keyword kw: Keyword arguments are passed on to the backend logging 
            system.
                        
        """
        kwargs = self._get_logging_kwargs(namespace, **kw)
        self.backend.critical(msg, **kwargs)

    def fatal(self, msg, namespace=None, **kw):
        """
        Log to the FATAL facility.
        
        :param msg: The message the log.
        :param namespace: A log prefix, generally the module ``__name__`` that 
            the log is coming from.  Will default to self._meta.namespace if 
            None is passed.  For debugging, it can be useful to set this to 
            ``__file__``, though ``__name__`` is much less verbose.
        :keyword kw: Keyword arguments are passed on to the backend logging 
            system.
                        
        """
        kwargs = self._get_logging_kwargs(namespace, **kw)
        self.backend.fatal(msg, **kwargs)

    def error(self, msg, namespace=None, **kw):
        """
        Log to the ERROR facility.
        
        :param msg: The message the log.
        :param namespace: A log prefix, generally the module ``__name__`` that 
            the log is coming from.  Will default to self._meta.namespace if 
            None is passed.  For debugging, it can be useful to set this to 
            ``__file__``, though ``__name__`` is much less verbose.
        :keyword kw: Keyword arguments are passed on to the backend logging 
            system.
                        
            """
        kwargs = self._get_logging_kwargs(namespace, **kw)
        self.backend.error(msg, **kwargs)
    
    ## NOTE: do we even need this for logbook?
    def clear_loggers(self):
        """Clear any previously configured logging namespaces.
        """
        if not self._meta.namespace:
            # _setup() probably wasn't run
            return
        self.handlers = []
