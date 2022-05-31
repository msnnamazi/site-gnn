class Error(Exception):
    """Base class for other exceptions"""
    pass


class PDBQTGenerationError(Error):
    """Raised when pdbqt generation failed."""
    def __init__(self):
        super().__init__("pdbqt generation failed.")


class PDBQTReadingError(Error):
    """Raised when pdbqt reading failed"""
    def __init__(self):
        super().__init__("pdbqt reading failed.")