import logging

# Clear all handlers from the root logger
logging.root.handlers.clear()

# Configure the logger
logger = logging.getLogger("vasp-input")
logger.setLevel(logging.INFO)
logger.propagate = False

# Ensure no duplicate handlers
if not logger.handlers:
    # Create a console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)

    # Define a formatter
    formatter = logging.Formatter("%(levelname)s %(message)s")
    console_handler.setFormatter(formatter)

    # Add the handler to the logger
    logger.addHandler(console_handler)

# Check if any other handlers are present
#print(f"Handlers attached to 'vasp-input': {logger.handlers}")
