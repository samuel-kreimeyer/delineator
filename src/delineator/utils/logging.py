"""Logging setup for the delineator package."""

import logging
import sys
from typing import Optional

# Package logger
_logger: Optional[logging.Logger] = None


def setup_logging(verbose: bool = False) -> None:
    """Set up logging for the package.

    Args:
        verbose: If True, set level to DEBUG; otherwise INFO.
    """
    global _logger

    level = logging.DEBUG if verbose else logging.INFO

    # Create formatter
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Create console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(level)
    console_handler.setFormatter(formatter)

    # Get package logger
    _logger = logging.getLogger("delineator")
    _logger.setLevel(level)

    # Remove existing handlers to avoid duplicates
    _logger.handlers.clear()
    _logger.addHandler(console_handler)


def get_logger(name: str) -> logging.Logger:
    """Get a logger for the specified module.

    Args:
        name: Module name (usually __name__).

    Returns:
        Logger instance.
    """
    # Ensure we have a child of the package logger
    if name.startswith("delineator"):
        return logging.getLogger(name)
    else:
        return logging.getLogger(f"delineator.{name}")
