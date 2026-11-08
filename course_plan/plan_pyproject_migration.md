# 計画：requirements.txt から pyproject.toml + uv.lock への移行

> **状態**：計画（未着手）．記録日 2026-11-08．
> 移行作業そのものは別途実施する．

## Context

現在の依存関係は `course/requirements.txt`（`uv pip freeze` の出力，142 パッケージ）で管理している．直接使うパッケージと推移的依存の区別がなく，追加・削除時に不要な依存が残る．Python の版は `course/.python-version` にしかない．今後モジュール 6 発展版で pertpy・PyDESeq2（推移的依存を含め 29 パッケージ）を追加する予定であり，その前に uv のプロジェクト管理（`pyproject.toml` + `uv.lock`）へ移行する．

決定事項（ユーザー回答）：

- `pyproject.toml` と `uv.lock` は **リポジトリ直下** に置く．`.venv` もリポジトリ直下に移す．
- `requirements.txt` は **削除** する．

## 方針

### 1. `pyproject.toml`（リポジトリ直下，新規）

直接使うパッケージだけを，現在の版を下限として書く．

```toml
[project]
name = "study-sc-rna-seq"
version = "0.1.0"
description = "数学者向け scRNA-seq 解析入門ハンズオン（Scanpy / PBMC 3k）"
requires-python = ">=3.13"
dependencies = [
    "scanpy>=1.12.1",
    "anndata>=0.12.11",
    "numpy>=2.4.4",
    "pandas>=2.3.3",
    "scipy>=1.17.1",
    "matplotlib>=3.10.8",
    "umap-learn>=0.5.12",
    "leidenalg>=0.11.0",
    "igraph>=1.0.0",
    "scvelo>=0.3.4",
    "jupyterlab>=4.5.6",
    "ipykernel>=7.2.0",
    "ipywidgets",          # tqdm の進捗表示（scVelo）を現状どおり保つ
    "nbconvert>=7.17.1",
]

[tool.uv]
package = false            # 配布用パッケージではない
```

選定根拠：ノートブックの import は `numpy, pandas, scanpy, scipy.sparse, scvelo`．加えて Leiden（`flavor='igraph'`）に `igraph`/`leidenalg`，UMAP に `umap-learn`，実行・閲覧に `jupyterlab`/`ipykernel`/`nbconvert`．`jupyter` メタパッケージと `python-igraph`（`igraph` の別名）は外す．

モジュール 6 発展版の pertpy・PyDESeq2 は，発展版の作業時に依存グループとして追加する（`uv add --group multidonor pertpy pydeseq2`）．本移行には含めない．

### 2. 現在の版を保ったままロックする

`>=` 下限だけで `uv lock` すると，最新版に上がりうる．移行で結果が変わらないよう，一時的に旧 `requirements.txt` を制約として与えてロックする：

1. `pyproject.toml` に一時的に `[tool.uv] constraint-dependencies = [...]`（旧 `requirements.txt` の全行）を入れて `uv lock`．
2. 制約を削除して再度 `uv lock`．uv は既存ロックの版を優先するので版は変わらない．
3. `uv export --format requirements-txt --no-hashes` と旧 `requirements.txt` を比べ，共通パッケージの版が一致することを確認する（`jupyter` メタパッケージ由来の一部パッケージは消えてよい）．

### 3. ファイルの移動・削除

| ファイル | 操作 |
|---|---|
| `pyproject.toml` | 新規（リポジトリ直下） |
| `uv.lock` | 新規（リポジトリ直下，コミットする） |
| `.python-version` | `course/` からリポジトリ直下へ移動（`3.13`） |
| `course/requirements.txt` | 削除 |
| `course/.venv` | 削除し，`uv sync` でリポジトリ直下に `.venv` を作る（`.gitignore` の `.venv` で除外済み） |

ノートブックのデータパスは `course/` 基準のまま変えない（Jupyter/nbconvert はカーネルをノートブックのあるディレクトリで起動するため）．

### 4. 文書の更新

- `README.md`：
  - 構成図：`requirements.txt` を削除し，直下の `pyproject.toml`・`uv.lock`・`.python-version` を追加．
  - 環境構築：`uv sync`（リポジトリ直下）のみに．
  - 実行：`uv run jupyter lab course/sc-rna-seq.ipynb`．
  - 主なパッケージの版の記述は `uv.lock` を参照する形にする．
  - データ取得手順（`cd course` してダウンロード）は変更なし．
- `CLAUDE.md`：
  - ノートブック実行：`uv run jupyter nbconvert --to notebook --execute --inplace --allow-errors --ExecutePreprocessor.timeout=1800 course/sc-rna-seq.ipynb`（リポジトリ直下から）．
  - パッケージ管理：「`uv pip freeze > requirements.txt`」を「`uv add <pkg>`（発展版用は `uv add --group <名前>`）で追加し，`pyproject.toml` と `uv.lock` をコミットする」に置き換える．
  - Python の版の記述を直下の `.python-version` に．
- `course_plan/scrna_hands_on_course.md`（モジュール 0 環境構築）：`uv venv` + `uv pip install ...` を `uv sync` に更新．
- `course_plan/plan_module6_multidonor.md`：「`requirements.txt` に反映」を「`uv add --group multidonor` で追加」に更新．

## 検証

1. `uv lock --check`（ロックが `pyproject.toml` と整合）．
2. 旧 `requirements.txt` と `uv export` の版比較（上記 2-3）．
3. `uv sync` 後，`uv run python -c "import scanpy, scvelo, leidenalg, igraph"` と版の表示．
4. `uv run jupyter nbconvert ... course/sc-rna-seq.ipynb` で全体を実行し，エラーがないこと，文字出力（HVG 数，発現割合表，擬似バルク表と検算）が移行前と一致することを確認（前回と同じ比較スクリプト）．
5. `uv run jupyter lab` が起動し，カーネルが直下の `.venv` を使うことを確認．
6. `git status` で，`requirements.txt` の削除・`.python-version` の移動・新規ファイルが意図どおりであることを確認．

## コミット

- 計画書：`course_plan/plan_pyproject_migration.md` として 1 コミット（今回）．
- 移行本体：上記をまとめて 1 コミット（別途，依頼時）．
