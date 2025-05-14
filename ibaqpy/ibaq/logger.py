import logging
import os
import sys
from typing import Optional, Dict, Any, Union
from datetime import datetime

# Default log format
DEFAULT_LOG_FORMAT = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
DEFAULT_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"

# Log levels dictionary for easy conversion from string
LOG_LEVELS = {
    "debug": logging.DEBUG,
    "info": logging.INFO,
    "warning": logging.WARNING,
    "error": logging.ERROR,
    "critical": logging.CRITICAL,
}


class ContextAdapter(logging.LoggerAdapter):
    """
    A logger adapter that adds context information to log messages.
    This makes logs more useful for debugging by providing additional context.
    """
    def process(self, msg, kwargs):
        if self.extra:
            context_str = " ".join(f"{k}={v}" for k, v in self.extra.items())
            return f"{msg} [{context_str}]", kwargs
        return msg, kwargs


def get_logger(name: str, context: Optional[Dict[str, Any]] = None) -> Union[logging.Logger, logging.LoggerAdapter]:
    """
    Get a logger with the specified name and optional context.
    
    Args:
        name: The name of the logger
        context: Optional dictionary of context values to include in log messages
        
    Returns:
        A logger or logger adapter with the specified name and context
    """
    logger = logging.getLogger(name)
    
    # If context is provided, return a ContextAdapter
    if context:
        return ContextAdapter(logger, context)
    
    return logger


def configure_logging(
    level: str = "info",
    log_file: Optional[str] = None,
    log_format: str = DEFAULT_LOG_FORMAT,
    date_format: str = DEFAULT_DATE_FORMAT,
    propagate: bool = True,
) -> None:
    """
    Configure the logging system for the application.
    
    Args:
        level: The log level (debug, info, warning, error, critical)
        log_file: Optional path to a log file
        log_format: The format string for log messages
        date_format: The format string for timestamps
        propagate: Whether to propagate logs to parent loggers
    """
    # Convert level string to logging level
    log_level = LOG_LEVELS.get(level.lower(), logging.INFO)
    
    # Configure root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(log_level)
    
    # Remove existing handlers to avoid duplicate logs
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)
    
    # Create formatter
    formatter = logging.Formatter(log_format, date_format)
    
    # Configure console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    root_logger.addHandler(console_handler)
    
    # Configure file handler if log_file is specified
    if log_file:
        os.makedirs(os.path.dirname(os.path.abspath(log_file)), exist_ok=True)
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        root_logger.addHandler(file_handler)
    
    # Configure ibaqpy loggers
    ibaqpy_logger = logging.getLogger("ibaqpy")
    ibaqpy_logger.setLevel(log_level)
    ibaqpy_logger.propagate = propagate


def log_execution_time(logger: Union[logging.Logger, logging.LoggerAdapter], level: int = logging.INFO):
    """
    Decorator to log the execution time of a function.
    
    Args:
        logger: The logger to use
        level: The log level to use
        
    Returns:
        A decorator that logs the execution time of the decorated function
    """
    def decorator(func):
        def wrapper(*args, **kwargs):
            start_time = datetime.now()
            logger.log(level, f"Starting {func.__name__}")
            
            try:
                result = func(*args, **kwargs)
                end_time = datetime.now()
                execution_time = end_time - start_time
                logger.log(level, f"Completed {func.__name__} in {execution_time}")
                return result
            except Exception as e:
                end_time = datetime.now()
                execution_time = end_time - start_time
                logger.exception(f"Error in {func.__name__} after {execution_time}: {str(e)}")
                raise
                
        return wrapper
    return decorator


def log_function_call(logger: Union[logging.Logger, logging.LoggerAdapter], level: int = logging.DEBUG):
    """
    Decorator to log function calls with arguments.
    
    Args:
        logger: The logger to use
        level: The log level to use
        
    Returns:
        A decorator that logs function calls with arguments
    """
    def decorator(func):
        def wrapper(*args, **kwargs):
            args_str = ", ".join([str(arg) for arg in args])
            kwargs_str = ", ".join([f"{k}={v}" for k, v in kwargs.items()])
            all_args = ", ".join(filter(None, [args_str, kwargs_str]))
            
            logger.log(level, f"Calling {func.__name__}({all_args})")
            return func(*args, **kwargs)
                
        return wrapper
    return decorator