from typing import Iterable


class SyntheticProfile():
    """
    simulated disease profile
    """

    def __init__(
            self,
            id: str,
            phenotypes: Iterable[str],
            disease: str = None):
        self.id = id
        self.phenotypes = phenotypes
        self.disease = disease
