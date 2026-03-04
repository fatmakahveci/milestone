from __future__ import annotations

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent


def test_static_site_assets_exist_and_contain_basic_seo_tags() -> None:
    index_html = (REPO_ROOT / "site" / "index.html").read_text(encoding="utf-8")
    benchmark_html = (REPO_ROOT / "site" / "benchmark.html").read_text(encoding="utf-8")
    robots = (REPO_ROOT / "site" / "robots.txt").read_text(encoding="utf-8")
    sitemap = (REPO_ROOT / "site" / "sitemap.xml").read_text(encoding="utf-8")

    assert "<title>Milestone | wgMLST web app and local workflow" in index_html
    assert 'name="description"' in index_html
    assert "benchmark packs provide reproducible wgMLST validation" in benchmark_html
    assert "Sitemap:" in robots
    assert "<urlset" in sitemap
