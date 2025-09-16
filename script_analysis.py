# -*- coding: utf-8 -*-
"""
script_analysis.py

Reproduz as análises do manuscrito "Impact of Digital Automation of Institutional 
Pathways on Test Ordering Patterns in a Primary Health Care Setting", conforme 
descrito no artigo (i.e., antes: 2023; depois: 2024). Trata-se de um script 
estruturado para replicabilidade e transparência, no padrão exigido pela revista 
Einstein (São Paulo).

SAÍDAS (na pasta --outdir):

  1) output_completo.xlsx
     - Contagens por exame em 2023 e 2024
     - Δ absoluto e Δ %
     - Flag "is_protocolo" (exame do protocolo institucional)
  2) output_completo_normalizado.xlsx
     - Mesmas colunas + "2023 per 100 visits", "2024 per 100 visits"
     - Δ absolute (per 100), Δ % (per 100)
  3) output_completo_normalizado_com_significancia.xlsx
     - Além das colunas acima: IRR (2024/2023), IC95% do IRR e p-valor
  4) summary_results.txt
     - Resumo dos principais números relatados no manuscrito
  5) (Opcional) slopegraph_top20.png
     - Gráfico 2023→2024 para os 20 exames mais frequentes em 2024

USO:

  python script_analysis.py \
      --data dataset_final_sample_20250908_1005.xlsx \
      --protocol de_para_nomes_exames_dataset.txt \
      --outdir ./outputs \
      --plot-slopegraph

Requisitos:

  Python versão 3.10+ (recomendado versao 3.11)
  Pacotes: pandas, numpy, scipy, openpyxl, matplotlib
"""

# ------------------------------------------------------------------------------
# Pacotes python
# ------------------------------------------------------------------------------

from __future__ import annotations

import argparse
import math
import unicodedata
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# SciPy apenas para o teste de Mann–Whitney (como no manuscrito)
from scipy.stats import mannwhitneyu

# ------------------------------------------------------------------------------
# Helpers gerais
# ------------------------------------------------------------------------------

def normalize_text(s: str) -> str:
    """Normaliza nomes de exames:
    - strip + upper
    - remove acentos (NFKD)
    - colapsa múltiplos espaços
    """

    if not isinstance(s, str):

        return s

    s = s.strip().upper()
    s = unicodedata.normalize("NFKD", s)
    s = "".join(ch for ch in s if not unicodedata.combining(ch))
    s = " ".join(s.split())

    return s


def norm_cdf(x: float) -> float:
    """Função de distribuição acumulada normal padrão."""

    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def irr_pvalue_ci(c1: float, e1: float, c2: float, e2: float) -> Tuple[float, Tuple[float, float], float, bool]:
    """Razão de taxas de Poisson (IRR), com correção de Haldane–Anscombe quando há zero.
    - IRR = (c2/e2) / (c1/e1)
    - SE(ln IRR) ≈ sqrt(1/c1 + 1/c2)  [aproximação clássica]
    - z = ln(IRR) / SE
    - p-valor bicaudal = 2 * (1 - Phi(|z|))
    - IC95% = exp( ln(IRR) ± 1.96*SE )

    Retorna: (irr, (lo, hi), p, zero_corrected_flag)
    """

    if (c1 + c2) == 0:
    
        return float('nan'), (float('nan'), float('nan')), float('nan'), False

    # Correção de Haldane–Anscombe se algum count==0
    c1a = c1 + (0.5 if c1 == 0 else 0.0)
    c2a = c2 + (0.5 if c2 == 0 else 0.0)

    irr = (c2a / e2) / (c1a / e1)
    se = math.sqrt((1.0 / c1a) + (1.0 / c2a))

    if se == 0 or not math.isfinite(irr) or irr <= 0:
    
        return float('nan'), (float('nan'), float('nan')), float('nan'), (c1 == 0 or c2 == 0)

    lnirr = math.log(irr)
    z = lnirr / se
    p = 2.0 * (1.0 - norm_cdf(abs(z)))
    lo = math.exp(lnirr - 1.96 * se)
    hi = math.exp(lnirr + 1.96 * se)
    zero_corrected = (c1 == 0 or c2 == 0)

    return irr, (lo, hi), p, zero_corrected


def shannon_entropy(probs: np.ndarray, base: float = 2.0) -> float:
    """Entropia de Shannon: H = -∑ p log_b p (ignorando p=0)."""

    p = probs[probs > 0]

    if p.size == 0:

        return float('nan')

    logs = np.log(p) / np.log(base)

    return float(-np.sum(p * logs))


def perplexity_from_entropy(H: float, base: float = 2.0) -> float:
    """Perplexidade = b^H (interpretação: número efetivo de categorias equiprováveis)."""

    if not math.isfinite(H):

        return float('nan')

    return float(base ** H)


@dataclass
class StudyPaths:

    data: Path
    protocol: Path
    outdir: Path


# ------------------------------------------------------------------------------
# Funções principais da análise
# ------------------------------------------------------------------------------

def load_inputs(paths: StudyPaths) -> Tuple[pd.DataFrame, pd.Series]:
    """Carrega dataset e a lista de exames do protocolo, aplicando normalização."""

    df = pd.read_excel(paths.data)

    with open(paths.protocol, "r", encoding="utf-8") as f:

        proto_list_raw = [ln.strip() for ln in f if ln.strip()]

    proto_unique_norm = sorted({normalize_text(x) for x in proto_list_raw})

    # Normalização no dataset
    df["exam_name_upper"] = df["exam_name"].apply(normalize_text)
    df["is_protocolo"] = df["exam_name_upper"].isin(set(proto_unique_norm))

    return df, pd.Series(proto_unique_norm, name="protocol_exam_names")


def build_tables(df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[int, int]]:
    """Gera:
      - tabela_completa_sorted (contagens por exame × ano + deltas)
      - tabela_norm (normalizada por 100 atendimentos + deltas)
      - dicionário com nº de atendimentos por ano
    """

    # Pivot básico
    base = (
        df.groupby(["year", "exam_name_upper", "is_protocolo"])
          .size()
          .reset_index(name="count")
          .pivot(index=["exam_name_upper", "is_protocolo"], columns="year", values="count")
          .fillna(0)
          .astype(int)
          .reset_index()
    )

    # Garante ambas as colunas de ano
    for y in (2023, 2024):

        if y not in base.columns:
    
            base[y] = 0

    # Deltas
    base["Δ absoluto"] = base[2024] - base[2023]
    den = base[2023].replace(0, 1)  # evita div/0 para Δ%
    base["Δ %"] = ((base[2024] - base[2023]) / den) * 100

    tabela_completa_sorted = base.sort_values(by=2024, ascending=False).reset_index(drop=True)

    # Exposição (nº de atendimentos únicos)
    atendimentos_por_ano = df.groupby("year")["encounter_id"].nunique().to_dict()
    n2023 = atendimentos_por_ano.get(2023, 1) or 1
    n2024 = atendimentos_por_ano.get(2024, 1) or 1

    # Normalizado por 100 atendimentos
    tabela_norm = tabela_completa_sorted.copy()
    tabela_norm["2023 per 100 visits"] = (tabela_norm[2023] / n2023) * 100
    tabela_norm["2024 per 100 visits"] = (tabela_norm[2024] / n2024) * 100

    den_tx = tabela_norm["2023 per 100 visits"].replace(0, 1)
    tabela_norm["Δ absolute (per 100)"] = tabela_norm["2024 per 100 visits"] - tabela_norm["2023 per 100 visits"]
    tabela_norm["Δ % (per 100)"] = ((tabela_norm["2024 per 100 visits"] - tabela_norm["2023 per 100 visits"]) / den_tx) * 100

    return tabela_completa_sorted, tabela_norm, atendimentos_por_ano


def add_significance(tabela: pd.DataFrame, exposure: Dict[int, int]) -> pd.DataFrame:
    """Adiciona IRR (2024/2023), IC95% e p-valor por exame."""

    e1, e2 = exposure.get(2023, 1) or 1, exposure.get(2024, 1) or 1

    irr_list, lo_list, hi_list, p_list, zero_corr = [], [], [], [], []

    for _, r in tabela.iterrows():

        c1, c2 = int(r[2023]), int(r[2024])
        irr, (lo, hi), p, zc = irr_pvalue_ci(c1, e1, c2, e2)
        irr_list.append(irr)
        lo_list.append(lo)
        hi_list.append(hi)
        p_list.append(p)
        zero_corr.append("Yes" if zc else "No")

    out = tabela.copy()
    out["IRR (2024/2023)"] = irr_list
    out["IRR 95% CI (lo)"] = lo_list
    out["IRR 95% CI (hi)"] = hi_list
    out["p-value (rate ratio)"] = p_list
    out["Signif. 5%"] = out["p-value (rate ratio)"].apply(
        lambda x: "Yes" if (pd.notnull(x) and x < 0.05) else ("NA" if pd.isnull(x) else "No")
    )
    out["Zero correction?"] = zero_corr

    return out


def per_encounter_adherence(df: pd.DataFrame) -> Dict[str, float | Tuple[float, float]]:
    """Calcula aderência por atendimento (contagem de exames do protocolo por encontro) e
    executa Mann–Whitney (two-sided). Retorna resumo com n, mediana, IQR, média±dp e p-valor."""

    mask = df["is_protocolo"]

    # Contagens por atendimento
    by_enc_year = (df[mask]
                   .groupby(["year", "encounter_id"])
                   .size()
                   .reset_index(name="adherent_count"))

    # Inclui zeros
    all_enc = df[["year", "encounter_id"]].drop_duplicates()
    per_enc_full = all_enc.merge(by_enc_year, on=["year", "encounter_id"], how="left")
    per_enc_full["adherent_count"] = per_enc_full["adherent_count"].fillna(0).astype(int)

    x2023 = per_enc_full.loc[per_enc_full["year"] == 2023, "adherent_count"].to_numpy()
    x2024 = per_enc_full.loc[per_enc_full["year"] == 2024, "adherent_count"].to_numpy()

    # Mann–Whitney
    mw = mannwhitneyu(x2023, x2024, alternative="two-sided")

    summary = {
        "n_2023": int(x2023.size),
        "n_2024": int(x2024.size),
        "median_2023": float(np.median(x2023)),
        "median_2024": float(np.median(x2024)),
        "IQR_2023": (float(np.percentile(x2023, 25)), float(np.percentile(x2023, 75))),
        "IQR_2024": (float(np.percentile(x2024, 25)), float(np.percentile(x2024, 75))),
        "mean_sd_2023": (float(np.mean(x2023)), float(np.std(x2023, ddof=1))),
        "mean_sd_2024": (float(np.mean(x2024)), float(np.std(x2024, ddof=1))),
        "U_statistic": float(mw.statistic),
        "p_value": float(mw.pvalue),
    }
    return summary


def entropy_perplexity(df: pd.DataFrame, base: float = 2.0) -> Dict[int, Dict[str, float]]:
    """Calcula entropia/perplexidade por ano com base na distribuição de contagens por exame."""

    dist = (df.groupby(["year", "exam_name_upper"])
              .size()
              .reset_index(name="count"))

    def _calc(year: int) -> Dict[str, float]:
    
        sub = dist[dist["year"] == year].copy()
        total = sub["count"].sum()
        probs = (sub["count"] / total).to_numpy()
        H = shannon_entropy(probs, base=base)
        P = perplexity_from_entropy(H, base=base)
    
        return {"H": H, "perplexity": P, "total": float(total), "distinct_tests": float(sub.shape[0])}

    return {y: _calc(y) for y in (2023, 2024)}


def top10_rankings(tabela: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Top-10 por Δ absoluto e Top-10 por Δ % (ordenação desc.)."""
    top_abs = (tabela
               .sort_values(by="Δ absoluto", ascending=False)
               .head(10)
               .reset_index(drop=True))
    top_pct = (tabela
               .sort_values(by="Δ %", ascending=False)
               .head(10)
               .reset_index(drop=True))
    return top_abs, top_pct


def save_outputs(
    outdir: Path,
    tabela_completa_sorted: pd.DataFrame,
    tabela_norm: pd.DataFrame,
    tabela_sig: pd.DataFrame,
    summary_txt: str,
) -> None:
    """Salva planilhas e resumo de resultados."""

    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "output_completo.xlsx").write_bytes(
        tabela_completa_sorted.to_excel(index=False).encode() if hasattr(pd.DataFrame, "to_excel_str") else b""
    )  # fallback protegido abaixo

    # Escreve de forma robusta (openpyxl)
    tabela_completa_sorted.to_excel(outdir / "output_completo.xlsx", index=False)
    tabela_norm.to_excel(outdir / "output_completo_normalizado.xlsx", index=False)
    tabela_sig.to_excel(outdir / "output_completo_normalizado_com_significancia.xlsx", index=False)

    with open(outdir / "summary_results.txt", "w", encoding="utf-8") as f:

        f.write(summary_txt)


def maybe_plot_slopegraph(tabela: pd.DataFrame, outdir: Path, do_plot: bool) -> None:
    """Slopegraph simples para top 20 por contagem em 2024."""

    if not do_plot:

        return

    try:

        plot_df = tabela.head(20).copy()

        plt.figure(figsize=(6, 10))

        for _, r in plot_df.iterrows():

            plt.plot([2023, 2024], [r[2023], r[2024]], marker="o", linewidth=1)
    
        plt.title("Exam counts per year (top 20 by 2024 count)")
        plt.xlabel("Year")
        plt.ylabel("Number of orders")
        plt.xticks([2023, 2024])
        plt.grid(axis="y", linestyle="--", alpha=0.3)
        plt.tight_layout()
        plt.savefig(outdir / "slopegraph_top20.png", dpi=200)
        plt.close()

    except Exception as exc:

        print(f"[warn] slopegraph não gerado: {exc}")


# ------------------------------------------------------------------------------
# Função principal (pipeline)
# ------------------------------------------------------------------------------

def main(paths: StudyPaths, make_plot: bool = False) -> None:

    # 1) Carrega dados e protocolo, normaliza nomes
    df, proto_series = load_inputs(paths)

    # 2) Reconstrói tabelas (contagens e normalizadas) + exposição
    tabela_completa_sorted, tabela_norm, atendimentos_por_ano = build_tables(df)

    # 3) Significância (IRR/IC95/p) por exame com base na exposição observada
    tabela_sig = add_significance(tabela_completa_sorted, atendimentos_por_ano)

    # 4) Aderência por atendimento (Mann–Whitney)
    per_enc_summary = per_encounter_adherence(df)

    # 5) Entropia e perplexidade por ano (número efetivo de exames)  # ver conceito no manuscrito  :contentReference[oaicite:2]{index=2}
    ep = entropy_perplexity(df, base=2.0)

    # 6) Top-10 (Δ abs e Δ %)
    top_abs, top_pct = top10_rankings(tabela_completa_sorted)

    # 7) Texto de resumo (espelha a apresentação do manuscrito)
    s = []
    s.append("=== Replication summary ===\n")
    s.append(f"Unique encounters per year: {atendimentos_por_ano}\n")
    s.append(f"Per-encounter adherence (Mann–Whitney p-value): {per_enc_summary['p_value']:.3f}\n")
    s.append(f"Medians (2023→2024): {per_enc_summary['median_2023']} → {per_enc_summary['median_2024']}  "
             f"IQR 2023={per_enc_summary['IQR_2023']} | IQR 2024={per_enc_summary['IQR_2024']}\n")
    s.append(f"Entropy/Perplexity:\n"
             f"  2023: H={ep[2023]['H']:.3f} bits, perplexity={ep[2023]['perplexity']:.1f} "
             f"(distinct={ep[2023]['distinct_tests']:.0f}, total={ep[2023]['total']:.0f})\n"
             f"  2024: H={ep[2024]['H']:.3f} bits, perplexity={ep[2024]['perplexity']:.1f} "
             f"(distinct={ep[2024]['distinct_tests']:.0f}, total={ep[2024]['total']:.0f})\n")
    s.append("\nTop-10 by absolute growth (Δ absoluto):\n")
    s.append(top_abs[["exam_name_upper", "is_protocolo", "Δ absoluto", 2023, 2024]].to_string(index=False))
    s.append("\n\nTop-10 by percentage growth (Δ %):\n")
    s.append(top_pct[["exam_name_upper", "is_protocolo", "Δ %", 2023, 2024]].to_string(index=False))
    summary_txt = "\n".join(s)

    # 8) Salva tudo
    paths.outdir.mkdir(parents=True, exist_ok=True)
    tabela_completa_sorted.to_excel(paths.outdir / "output_completo.xlsx", index=False)
    tabela_norm.to_excel(paths.outdir / "output_completo_normalizado.xlsx", index=False)
    tabela_sig.to_excel(paths.outdir / "output_completo_normalizado_com_significancia.xlsx", index=False)

    with open(paths.outdir / "summary_results.txt", "w", encoding="utf-8") as f:

        f.write(summary_txt)

    # 9) (Opcional) Slopegraph
    maybe_plot_slopegraph(tabela_completa_sorted, paths.outdir, make_plot)

    # 10) Logs finais de sucesso
    print("✓ Outputs gravados em:", paths.outdir.resolve())
    print("  - output_completo.xlsx")
    print("  - output_completo_normalizado.xlsx")
    print("  - output_completo_normalizado_com_significancia.xlsx")
    print("  - summary_results.txt")
    if make_plot:

        print("  - slopegraph_top20.png")


# ------------------------------------------------------------------------------
# INTERFACE
# ------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Script para replicar as análises do manuscrito (einstein, SSDC/pathways)."
    )
    parser.add_argument("--data", type=Path, required=True,
                        help="Arquivo de dados (XLSX) — e.g., dataset_final_sample_20250908_1005.xlsx")
    parser.add_argument("--protocol", type=Path, required=True,
                        help="Arquivo TXT com os nomes de exames do protocolo (um por linha).")
    parser.add_argument("--outdir", type=Path, required=True,
                        help="Pasta de saída para gerar planilhas/relatórios.")
    parser.add_argument("--plot-slopegraph", action="store_true",
                        help="Se indicado, gera um slopegraph (top 20 por contagem em 2024).")

    args = parser.parse_args()
    paths = StudyPaths(data=args.data, protocol=args.protocol, outdir=args.outdir)

    # Verificações básicas
    if not paths.data.exists():

        raise FileNotFoundError(f"Arquivo de dados não encontrado: {paths.data}")

    if not paths.protocol.exists():

        raise FileNotFoundError(f"Arquivo de protocolo não encontrado: {paths.protocol}")

    main(paths, make_plot=args.plot_slopegraph)
