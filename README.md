# Impacto da Automação Digital de Protocolos Institucionais na Solicitação de Exames

Este repositório contém o código e os dados utilizados para as análises do manuscrito:

**"Impact of Digital Automation of Institutional Pathways on Test Ordering Patterns in a Primary Health Care Setting"**

O objetivo é avaliar como a automação digital de protocolos institucionais influencia o padrão de solicitação de exames em um serviço de Atenção Primária à Saúde.

## Estrutura do Projeto

- `script_analysis.py`: Script principal para análise dos dados e geração dos resultados.
- `dataset_final_sample_20250908_1005.xlsx`: Base de dados utilizada na análise.
- `de_para_nomes_exames_dataset.txt`: Mapeamento de nomes de exames para padronização.
- `outputs/`: Pasta com todos os resultados gerados pelo script.
	- `output_completo.xlsx`: Contagens por exame em 2023 e 2024.
	- `output_completo_normalizado.xlsx`: Contagens normalizadas por 100 atendimentos.
	- `output_completo_normalizado_com_significancia.xlsx`: Resultados com significância estatística.
	- `summary_results.txt`: Resumo dos principais resultados.
	- `slopegraph_top20.png`: Gráfico dos 20 exames mais frequentes.

## Como Executar

Recomenda-se utilizar Python 3.10+ (preferencialmente 3.11). Instale os pacotes necessários:

```bash
pip install pandas numpy scipy openpyxl matplotlib
```

Execute o script principal:

```bash
python script_analysis.py \
		--data dataset_final_sample_20250908_1005.xlsx \
		--protocol de_para_nomes_exames_dataset.txt \
		--outdir ./outputs \
		--plot-slopegraph
```

## Resultados

Os resultados serão salvos na pasta `outputs/` conforme descrito acima.

## Contato

Para dúvidas ou sugestões, entre em contato com os autores do projeto.