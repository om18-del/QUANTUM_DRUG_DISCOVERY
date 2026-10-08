from __future__ import annotations

from app.config import Settings, project_root


def test_project_root_is_repo_root():
    settings = Settings()
    repo_root = project_root()
    settings_paths = [
        settings.cache_dir,
        settings.job_store_path,
        settings.h2_hamiltonian_dir,
        settings.complexes_manifest_path,
        settings.models_manifest_path,
    ]
    for path in settings_paths:
        resolved = path.resolve()
        assert repo_root in resolved.parents or resolved == repo_root, (
            f"Path {path} is outside repo root {repo_root}"
        )


def test_env_file_points_at_project_root():
    settings = Settings()
    expected_env = project_root() / ".env"
    assert settings.model_config["env_file"] == str(expected_env)
