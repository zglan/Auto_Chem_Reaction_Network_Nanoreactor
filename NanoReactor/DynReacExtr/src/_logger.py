import os
import sys
import logging


PATH = os.path.abspath('.') + '/logs/'
FMT = '%(asctime)s %(filename)s  %(levelname)s: %(message)s'
DATEFMT = '%Y-%m-%d %H:%M:%S'

class MyLog(object):
    def __init__(self):
        self.logger = logging.getLogger()
        self.formatter = logging.Formatter(fmt=FMT, datefmt=DATEFMT)

        self.logger.addHandler(self.get_console_handler())
        self.logger.setLevel(logging.INFO)
        self.logger.handlers.pop()

    def get_file_handler(self, filename):
        filehandler = logging.FileHandler(filename, encoding="utf-8")
        filehandler.setFormatter(self.formatter)
        return filehandler

    def get_console_handler(self):
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setFormatter(self.formatter)
        return console_handler
