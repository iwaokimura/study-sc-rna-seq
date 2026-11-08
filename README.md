# study-sc-rna-seq

数学者向けの scRNA-seq（single-cell RNA sequencing）解析入門ハンズオンの学習リポジトリ．

数理的素養と計算機操作の経験はあるが，分子生物学・医学生理学の背景知識をほとんど持たない人を対象に，定番ツール **Scanpy（Python）** と公開データ **10x Genomics PBMC 3k**（健常ドナー末梢血単核球，2,700 細胞）を用いて，典型的な解析パイプライン全体を体験する．各ステップでは「数理的に何をしているか」（行列演算，確率モデル，グラフ理論，多様体学習など）を重視する．

## 構成

```
study-sc-rna-seq/
├── pyproject.toml                 # 依存パッケージ（直接使うもの）
├── uv.lock                        # 依存関係のロック（全パッケージの版）
├── .python-version                # Python の版（3.13）
├── course_plan/
│   ├── sc-RNAseq-analysis.md      # 初版のコース設計
│   ├── scrna_hands_on_course.md   # 詳細版コースノート（最新）
│   └── plan_*.md                  # 今後の作業計画
└── course/
    ├── sc-rna-seq.ipynb           # ハンズオン実行ノートブック（モジュール 1〜7）
    ├── module6-multidonor.ipynb   # モジュール 6 発展版（複数ドナーの擬似バルク）
    └── data/                      # 解析データ（git 管理外．下記手順で取得）
```

## コース内容

| モジュール | 内容 | 対応する数理概念 |
|---|---|---|
| 0 | 準備：生物学の最小知識と環境構築 | — |
| 1 | カウント行列と品質管理（QC） | 非負整数値疎行列，外れ値検出 |
| 2 | 正規化・高変動遺伝子（HVG）選択 | スケーリング・対数変換，特徴選択 |
| 3 | 次元削減：PCA・UMAP | SVD，多様体学習 |
| 4 | クラスタリング（Leiden） | グラフのモジュラリティ最適化 |
| 5 | 細胞型アノテーション（マーカー遺伝子） | ノンパラメトリック検定（Wilcoxon） |
| 6 | 差次発現解析（DEG）・疑似バルク | 負の二項分布モデル |
| 7 | 発展：PAGA・拡散擬似時間・RNA velocity | マルコフ連鎖・拡散過程，常微分方程式 |

詳細な解説と参考文献は [`course_plan/scrna_hands_on_course.md`](course_plan/scrna_hands_on_course.md) を参照．

## 環境構築

[uv](https://docs.astral.sh/uv/) のプロジェクト管理（`pyproject.toml` + `uv.lock`）を用いる．Python は 3.13（`.python-version` で指定）．リポジトリ直下で

```bash
uv sync
```

を実行すると，`uv.lock` どおりの環境が `.venv/` に作られる．

主なパッケージ：scanpy 1.12.1，anndata 0.12.11，leidenalg 0.11.0，umap-learn 0.5.12，scvelo 0.3.4，jupyterlab 4.5.6（正確な版は `uv.lock` を参照）．

## データ取得

ノートブックは `course/` からの相対パス `data/filtered_gene_bc_matrices/hg19/` を読み込む．

```bash
cd course
mkdir -p data
curl https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz \
  -o data/pbmc3k_filtered_gene_bc_matrices.tar.gz
tar -xzf data/pbmc3k_filtered_gene_bc_matrices.tar.gz -C data/
```

モジュール 6 発展版では Kang et al. (2018) のデータを使う（`course/` で実行）：

```bash
curl -L https://exampledata.scverse.org/pertpy/kang_2018.h5ad -o data/kang_2018.h5ad
```

モジュール 7 の RNA velocity では，scVelo 付属のサンプルデータ（`scv.datasets.pancreas()`，膵臓内分泌細胞の分化）を実行時に自動ダウンロードする．

## 実行

リポジトリ直下で

```bash
uv run jupyter lab course/sc-rna-seq.ipynb
```

モジュール 6 発展版は PyDESeq2（依存グループ `multidonor`，既定の `uv sync` では入らない）を使うので，グループを指定して起動する．全体の実行は約 1 分である．

```bash
uv run --group multidonor jupyter lab course/module6-multidonor.ipynb
```

## 注意点

- PBMC 3k は単一ドナー由来のため，モジュール 6 の疑似バルク解析は，細胞を無作為に分けた仮サンプルでの集約の実習にとどめ，検定は行わない．複数ドナーのデータでの検定は発展版（`course/module6-multidonor.ipynb`）で行う．
- モジュール 7 の DPT で使う根細胞（`iroot`）は便宜的に指定したもの．PBMC は分化系列のデータではないため，擬似時間の生物学的解釈には注意が必要．

## 主な参考資料

- Scanpy PBMC 3k チュートリアル（legacy workflow）：https://scanpy.readthedocs.io/en/1.12.x/tutorials/basics/clustering-2017.html
- Single-cell best practices（Theis Lab）：https://www.sc-best-practices.org
- Luecken & Theis (2019), *Mol. Syst. Biol.* 15, e8746：https://doi.org/10.15252/msb.20188746
- Heumos et al. (2023), *Nat. Rev. Genet.* 24, 550–572：https://doi.org/10.1038/s41576-023-00586-w
