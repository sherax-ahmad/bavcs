from .score import (
    VariantScores,
    BAVCSResult,
    compute_bavcs,
    representation_confidence,
    VEP_CLASSES,
    POPULATION_FREE,
    POPULATION_TUNED,
    CLINICAL_TRAINED,
)
from .real_data import (
    load_real_sample,
    normalize_real_sample,
    compute_real_discordance,
    DiscordanceResult,
)

__all__ = [
    "VariantScores",
    "BAVCSResult",
    "compute_bavcs",
    "representation_confidence",
    "VEP_CLASSES",
    "POPULATION_FREE",
    "POPULATION_TUNED",
    "CLINICAL_TRAINED",
    "load_real_sample",
    "normalize_real_sample",
    "compute_real_discordance",
    "DiscordanceResult",
]

__version__ = "0.1.0"
