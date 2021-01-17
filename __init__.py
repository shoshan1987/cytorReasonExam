import os
import sys
import logging

from .commom.logging import set_logger

def init_logger(log_name, log_mode = 'a'):
    logger = logging.getLogger('logger')
    set_logger(logger, log_name, log_mode = log_mode)
    logger.info("started running")
    return(logger)