import os
from pathlib import Path

from tdpy.paths import RepositoryPaths


PATH_ENV_VAR = "AEGIS_PATH"
os.environ.setdefault(PATH_ENV_VAR, str(Path(__file__).resolve().parent))
_REPOSITORY_PATHS = RepositoryPaths(PATH_ENV_VAR)

get_repository_path = _REPOSITORY_PATHS.get_repository_path
get_data_path = _REPOSITORY_PATHS.get_data_path
get_visuals_path = _REPOSITORY_PATHS.get_visuals_path