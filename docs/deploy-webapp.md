# Publish Milestone Web Demo

The UI in `webapp/app.py` is a Streamlit app. It cannot run on GitHub Pages (Pages is static-only).

Use the production stack under `deploy/` to publish at:

- `https://<your-domain>/milestone/`

## 1. Prerequisites

- A Linux host (VM/VPS) with Docker + Docker Compose
- DNS for your domain pointing to that host's public IP
- Ports `80` and `443` open on the host firewall

## 2. Start the app

From the repository root:

```bash
cd deploy
export DOMAIN=fatmakahveci.com
docker compose -f docker-compose.prod.yml up -d --build
```

This stack runs:

- `milestone-web`: Streamlit app with `--server.baseUrlPath=milestone`
- `milestone-proxy`: Caddy reverse proxy + automatic HTTPS

## 3. Verify

```bash
docker compose -f docker-compose.prod.yml ps
docker compose -f docker-compose.prod.yml logs --tail=100 milestone-web
```

Open:

- `https://fatmakahveci.com/milestone/`

## 4. Updates

When you change code:

```bash
cd deploy
docker compose -f docker-compose.prod.yml up -d --build
```

## Notes

- Keep the GitHub Pages static landing page if you want, but it should link to the live app URL.
- If you already run another reverse proxy on the same server, do not run Caddy here; proxy `/milestone` to `milestone-web:8501` and keep Streamlit `baseUrlPath=milestone`.
