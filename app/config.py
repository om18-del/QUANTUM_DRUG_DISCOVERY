from __future__ import annotations

from functools import lru_cache
from pathlib import Path

from pydantic_settings import BaseSettings, SettingsConfigDict


def project_root() -> Path:
    return Path(__file__).resolve().parents[1]


class Settings(BaseSettings):
    model_config = SettingsConfigDict(
        env_file=str(project_root() / ".env"),
        env_file_encoding="utf-8",
        case_sensitive=False,
        extra="ignore",
    )

    # App
    app_name: str = "Quantum Drug Discovery Explorer API"
    api_prefix: str = "/api"

    # Server
    host: str = "127.0.0.1"
    port: int = 8000
    workers: int = 1

    # CORS
    cors_origins: str = ""
    cors_allow_credentials: bool = False

    # Cache
    cache_dir: Path = project_root() / "app" / "data" / "cache"
    cache_ttl_pubchem_hours: int = 24
    cache_ttl_pdb_hours: int = 24
    cache_ttl_smiles_descriptors_hours: int = 720

    # Job store
    job_store_path: Path = project_root() / "app" / "data" / "jobs.sqlite"

    # External API timeouts
    external_timeout_seconds: float = 10.0
    external_max_retries: int = 3

    # Input limits
    max_smiles_length: int = 200
    max_query_length: int = 200

    # Quantum
    h2_hamiltonian_dir: Path = project_root() / "app" / "data" / "hamiltonians"

    # Complexes
    complexes_manifest_path: Path = (
        project_root() / "app" / "data" / "curated_complexes" / "manifest.json"
    )

    # AI models
    models_manifest_path: Path = project_root() / "app" / "data" / "models" / "manifest.json"

    @property
    def cors_origin_list(self) -> list[str]:
        if not self.cors_origins.strip():
            return []
        return [o.strip() for o in self.cors_origins.split(",") if o.strip()]


@lru_cache
def get_settings() -> Settings:
    return Settings()
