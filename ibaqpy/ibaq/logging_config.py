"""
Logging configuration for the ibaqpy package.

This module provides functions to configure the logging system for the ibaqpy package.
It should be imported and initialized at the start of the application.
"""

import os
import logging
from typing import Optional

from ibaqpy.ibaq.logger import configure_logging


def initialize_logging(
    level: str = "info",
    log_file: Optional[str] = None,
    log_format: Optional[str] = None,
    date_format: Optional[str] = None,
) -> None:
    """
    Initialize the logging system for the ibaqpy package.

    This function should be called at the start of the application to configure
    the logging system. It sets up console and file logging with appropriate
    formatting.

    Args:
        level: The log level (debug, info, warning, error, critical)
        log_file: Optional path to a log file
        log_format: Optional format string for log messages
        date_format: Optional format string for timestamps
    """
    # Use environment variables if available
    env_level = os.environ.get("IBAQPY_LOG_LEVEL", level)
    env_log_file = os.environ.get("IBAQPY_LOG_FILE", log_file)

    # Configure logging
    configure_logging(
        level=env_level,
        log_file=env_log_file,
        log_format=log_format,
        date_format=date_format,
    )

    # Log initialization
    logger = logging.getLogger("ibaqpy")
    logger.info("Logging initialized at level %s", env_level.upper())
    if env_log_file:
        logger.info("Log file: %s", env_log_file)


def get_log_file_path(base_dir: Optional[str] = None) -> str:
    """
    Get a default log file path based on the current date.

    Args:
        base_dir: Optional base directory for log files

    Returns:
        A path to a log file
    """
    import datetime

    # Default to logs directory in current working directory
    if base_dir is None:
        base_dir = os.path.join(os.getcwd(), "logs")

    # Create logs directory if it doesn't exist
    os.makedirs(base_dir, exist_ok=True)

    # Create log file name based on current date
    date_str = datetime.datetime.now().strftime("%Y-%m-%d")
    log_file = os.path.join(base_dir, f"ibaqpy_{date_str}.log")

    return log_file