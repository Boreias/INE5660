#!/usr/bin/env python3
import os
import re
import csv

# Caminho base do projeto (ajuste se necessário)
BASE_DIR = "..\\Data\\Output\\"
OUTPUT_CSV = "dados_perf_consolidados.csv"

# Regex para capturar métricas
METRICS = {
    "cpu_clock": r"([\d\.,]+)\s+msec cpu-clock",
    "cycles": r"([\d\.,]+)\s+cycles",
    "cache_misses": r"([\d\.,]+)\s+cache-misses",
    "cache_refs": r"([\d\.,]+)\s+cache-references",
    "branch_misses": r"([\d\.,]+)\s+branch-misses",
    "branches": r"([\d\.,]+)\s+branches",
    "page_faults": r"([\d\.,]+)\s+page-faults",
    "time_elapsed": r"([\d\.,]+)\s+seconds time elapsed",
    "user_time": r"([\d\.,]+)\s+seconds user",
    "sys_time": r"([\d\.,]+)\s+seconds sys",
}

def parse_perf_file(filepath):
    """Extrai as métricas de um arquivo .txt do perf."""
    data = {}
    with open(filepath, "r", encoding="utf-8") as f:
        text = f.read()

    for key, pattern in METRICS.items():
        match = re.search(pattern, text)
        if match:
            # Substitui vírgulas por pontos e remove separadores de milhar
            value = match.group(1).replace(".", "").replace(",", ".")
            data[key] = float(value)
        else:
            data[key] = None

    return data

def extract_metadata_from_path(filepath):
    """Extrai linguagem, tipo, threads, número de entrada e caso do nome do arquivo."""
    parts = filepath.split(os.sep)
    info = {
        "linguagem": None,
        "tipo": None,
        "threads": None,
        "arquivo_entrada": None,
        "caso": None,
        "versao_do_codigo": None
    }

    # Detecta linguagem
    if "Rust" in parts:
        info["linguagem"] = "Rust"
    elif "Fortran" in parts:
        info["linguagem"] = "Fortran"

    # Detecta tipo
    if "Parallel" in parts:
        info["tipo"] = "Paralelo"
    else:
        info["tipo"] = "Sequencial"
    
    if "Duration" in parts:
        info["versao_do_codigo"] = "Duration"
    elif "General_without_Print" in parts:
        info["versao_do_codigo"] = "General_without_Print"
    else:
        info["versao_do_codigo"] = "General"

    filename = os.path.basename(filepath).replace(".txt", "")

    # --- Caso paralelo: THREAD_NUMERO_CASE.txt ---
    if info["tipo"] == "Paralelo":
        match = re.match(r"(\d+)_([\d_]+)_(\d+)$", filename)
        if match:
            info["threads"] = int(match.group(1))
            info["arquivo_entrada"] = match.group(2)
            info["caso"] = int(match.group(3))
        else:
            print(f"[!] Nome inesperado (paralelo): {filename}")

    # --- Caso sequencial: NUMERO_CASE.txt ---
    else:
        match = re.match(r"([\d_]+)_(\d+)$", filename)
        if match:
            info["arquivo_entrada"] = match.group(1)
            info["caso"] = int(match.group(2))
            info["threads"] = 1  # padrão: 1 thread
        else:
            print(f"[!] Nome inesperado (sequencial): {filename}")

    return info


def main():
    rows = []

    for root, _, files in os.walk(BASE_DIR):
        if root.endswith("Execution_Result_File"):
            for file in files:
                if file.endswith(".txt"):
                    path = os.path.join(root, file)
                    metrics = parse_perf_file(path)
                    if metrics.get('cpu_clock') != None:
                        meta = extract_metadata_from_path(path)
                        row = {**meta, **metrics, "arquivo": path}
                        rows.append(row)

    # Escreve o CSV
    with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as csvfile:
        fieldnames = [
            "linguagem", "tipo", "threads", "arquivo_entrada", "arquivo", "caso", "versao_do_codigo"
        ] + list(METRICS.keys())
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"[✔] Arquivo CSV gerado: {OUTPUT_CSV}")
    print(f"[ℹ] Total de registros: {len(rows)}")

if __name__ == "__main__":
    main()
