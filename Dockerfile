FROM python:3.11-slim

WORKDIR /app

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
ENV STREAMLIT_SERVER_PORT=8501
ENV STREAMLIT_SERVER_ADDRESS=0.0.0.0

COPY webapp/requirements.txt /tmp/webapp-requirements.txt
RUN pip install --no-cache-dir -r /tmp/webapp-requirements.txt

COPY . /app

EXPOSE 8501

CMD ["streamlit", "run", "webapp/app.py"]
