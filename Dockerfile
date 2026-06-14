FROM python:3.11-slim

WORKDIR /app

# Install build tools needed by scipy
RUN apt-get update && apt-get install -y --no-install-recommends \
        gcc \
        g++ \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies first (layer cache)
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application files
COPY scenario04-Cesium-advanced02.py .
COPY sat_metadata.csv .
COPY data/ ./data/
COPY space_db_slim.duckdb .

# HuggingFace Spaces requires port 7860
ENV PORT=7860
ENV HOST=0.0.0.0
ENV DB_PATH=/app/space_db_slim.duckdb

EXPOSE 7860

CMD ["python", "scenario04-Cesium-advanced02.py"]
