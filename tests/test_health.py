from __future__ import annotations

from fastapi.testclient import TestClient

from app.main import app


def test_health_returns_ok():
    with TestClient(app) as client:
        response = client.get("/health")
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "ok"
        assert "app" in data
        assert data["app"] == "Quantum Drug Discovery Explorer API"


def test_generic_error_returns_error_shape():
    with TestClient(app) as client:
        response = client.get("/nonexistent-route-for-error-test")
        assert response.status_code == 404
        data = response.json()
        assert "error" in data
        assert "code" in data["error"]
        assert "message" in data["error"]
