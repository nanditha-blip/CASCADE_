from flask import Flask, request, jsonify
from flask_cors import CORS
import os

from algorithms import analyze
from data import GENE_PATHWAYS

app = Flask(__name__)
CORS(app)  # allows the frontend (different origin / file) to call this API


@app.get("/api/health")
def health():
    return jsonify({"status": "ok"})


@app.get("/api/known-genes")
def known_genes():
    """Lets the frontend show an autocomplete / example list."""
    return jsonify({"genes": sorted(GENE_PATHWAYS.keys())})


@app.post("/api/analyze")
def analyze_genes():
    body = request.get_json(silent=True) or {}
    genes = body.get("genes", [])

    if not isinstance(genes, list) or not genes:
        return jsonify({"error": "Please provide a non-empty list of genes."}), 400
    if len(genes) > 100:
        return jsonify({"error": "Please limit uploads to 100 genes for this demo."}), 400

    result = analyze(genes)
    return jsonify(result)


if __name__ == "__main__":
    # host="0.0.0.0" so phones on the same wifi (or, once deployed, the
    # public internet) can reach this server. PORT is read from the
    # environment so hosting platforms like Render can assign it.
    port = int(os.environ.get("PORT", 5000))
    app.run(host="0.0.0.0", port=port, debug=True)
