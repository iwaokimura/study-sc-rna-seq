# scRNA-seq ハンズオンコース設計：数学者向け入門

## 概要

本コースは，数理的素養と計算機操作の経験を持つが，分子生物学・医学生理学の背景知識をほとんど持たない数学者を対象として設計した scRNA-seq（single-cell RNA sequencing）解析のハンズオンコースである．典型的な解析パイプライン全体を，定番ツール（Scanpy / Python）と公開データ（10x Genomics PBMC 3k データセット）を用いて体験することを目標とする．

生物学的概念は「数理的に何をしているか」という視点から必要最小限を導入し，各ステップの数学的本質（行列演算，確率モデル，グラフ理論，多様体学習など）を強調する．コース全体の所要時間の目安は約 **3〜5日間**（集中実施の場合）または隔週セッション×5〜6回である．

***

## コースの全体構成

| モジュール | タイトル | 所要時間目安 |
|---|---|---|
| 0 | 準備：生物学の最小知識と環境構築 | 0.5日 |
| 1 | データの正体を知る：カウント行列と QC | 0.5日 |
| 2 | 正規化・高変動遺伝子選択：統計モデルの話 | 0.5日 |
| 3 | 次元削減：PCA・UMAP の数理 | 0.5日 |
| 4 | クラスタリングとグラフ分割 | 0.5日 |
| 5 | 細胞型アノテーション：マーカー遺伝子と検定 | 0.5日 |
| 6 | 差次発現解析：負の二項分布モデル | 0.5日 |
| 7 | 発展：軌跡推定と擬似時間 | 0.5日 |

***

## モジュール 0：準備

### 最小限の生物学知識

数学者が scRNA-seq を理解するために必要な生物学的事実は以下の 3 点に絞られる：

1. **細胞** は生命の基本単位．ヒトの体には約 37 兆個の細胞があり，それぞれが同一のゲノム DNA を持つ．
2. **遺伝子発現**（転写）とは，DNA の情報が mRNA に写し取られるプロセス．mRNA の量が「その遺伝子がどれだけ使われているか」を示す．scRNA-seq はこの mRNA 量を各細胞ごとに計測する技術である．[^1]
3. **PBMC**（末梢血単核球）とは血液中の免疫細胞の集合体で，T 細胞（CD4+: 25〜60%，CD8+: 5〜30%），B 細胞（5〜10%，CD19+），NK 細胞（10〜30%，CD56+CD3-），単球（CD14+）などを含む．本コースで用いる PBMC データは，これらの異なる細胞種が混在する状況での分類問題として捉えると数学的に理解しやすい．[^2][^3]

> **数理的直観**：scRNA-seq の出力は \(m \times n\) 行列（\(m\) = 遺伝子数 ≈ 20,000，\(n\) = 細胞数 ≈ 3,000）であり，各エントリは非負整数（UMI カウント）である．目標はこの高次元空間の点群を意味のある部分集合（細胞型）に分類することである．

### 使用ツールと公開データ

- **メインツール**：[Scanpy](https://scanpy.readthedocs.io/) (Python)．AnnData オブジェクトをベースにした scRNA-seq 解析パッケージ[^4][^5]
- **データ**：10x Genomics の「PBMC 3k」（健常ドナー由来，2,700 細胞）．以下のコマンドで取得：[^4]

```bash
mkdir -p data write
curl https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz \
  -o data/pbmc3k_filtered_gene_bc_matrices.tar.gz
tar -xzf data/pbmc3k_filtered_gene_bc_matrices.tar.gz -C data/
```

- **参考書**：[Single-cell best practices](https://www.sc-best-practices.org/)（Theis ラボ，Helmholtz Munich）— 無料のオンライン教科書で，独立ベンチマーク研究に基づくベストプラクティスを随時更新[^6][^7]

### 環境構築

```bash
pip install scanpy anndata numpy pandas matplotlib seaborn umap-learn leidenalg python-igraph
```

Jupyter Notebook または Google Colab での作業を推奨する．[^8]

***

## モジュール 1：データの正体を知る — カウント行列と品質管理（QC）

### AnnData オブジェクト

```python
import scanpy as sc

adata = sc.read_10x_mtx('data/filtered_gene_bc_matrices/hg19/', var_names='gene_symbols')
adata  # AnnData object: n_obs(cells) × n_vars(genes)
```

`adata.X` は疎行列（scipy sparse matrix）であり，\(2{,}700 \times 32{,}738\) の UMI カウント行列である．エントリの大部分はゼロ（スパース性 > 90%）で，これは「ドロップアウト」と呼ばれる技術的アーティファクトに起因する．[^9]

### QC 指標の計算

```python
sc.pp.calculate_qc_metrics(adata, inplace=True)
# ミトコンドリア遺伝子の割合を追加
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
```

各細胞に対して計算する主要な QC 指標：

| 指標 | 数理的意味 | フィルタ閾値の目安 |
|---|---|---|
| `n_genes_by_counts` | 非ゼロエントリ数（\(\ell_0\) ノルム的） | 200〜2,500 |
| `total_counts` | 各行の \(\ell_1\) ノルム | 外れ値除去 |
| `pct_counts_mt` | ミトコンドリア遺伝子カウントの割合 | < 5% |

ミトコンドリア遺伝子の高発現は細胞の「死」の指標，遺伝子数が異常に多い細胞はダブレット（2 細胞が 1 つのバーコードに混入）の疑いがある．[^10][^9]

```python
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]
sc.pp.filter_genes(adata, min_cells=3)
```

> **発展**：アンビエント RNA（空の液滴に含まれる細胞外 mRNA の混入）除去には **SoupX** が有効である．ダブレット検出には **scDblFinder** が推奨される．本コースでは基本フィルタのみ実施する．[^11][^10]

***

## モジュール 2：正規化と高変動遺伝子選択

### 正規化：なぜ必要か

各細胞間で捕捉された総 RNA 量（ライブラリサイズ，\(s_i = \sum_j x_{ij}\)）は技術的ばらつきにより異なる．単純なサイズファクター正規化（CPM 正規化相当）を適用する：

\[ \tilde{x}_{ij} = \frac{x_{ij}}{s_i} \cdot 10{,}000 \]

その後対数変換 \(\log_1(1 + \tilde{x}_{ij})\) を施す（`log1p` 変換）．[^12][^7]

```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata  # 正規化済みデータを保存
```

より高度な正規化として **SCTransform**（Seurat）や Scran 正規化がある．SCTransform は細胞ごとのシーケンス深度を共変量としてピアソン残差を計算し，技術的ばらつきを除去する．[^12]

### 高変動遺伝子（HVG）の選択

全 32,000 遺伝子のうち，細胞間で統計的に有意な発現変動を示す遺伝子（高変動遺伝子, HVG）を選択して以降の解析に用いる．

```python
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
# 通常 2,000〜4,000 遺伝子が選択される
```

> **数理的注意点**：`sc.pp.scale()` で各遺伝子を標準化（平均 0，分散 1）することが多い．これは主成分分析の前処理として標準的だが，scRNA-seq では発現量 0 が多いため，スパース性を損なうというトレードオフがある．[^13]

***

## モジュール 3：次元削減 — PCA と UMAP の数理

### PCA：線形次元削減

正規化済みデータ行列 \(X \in \mathbb{R}^{n \times m}\)（\(n\) = 細胞，\(m\) = HVG 数）に対し，PCA は共分散行列の固有値分解 \(X^\top X = V \Lambda V^\top\) を通じて，分散最大方向への射影を求める．[^14]

```python
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True)  # elbow plot
```

通常上位 10〜50 PC を使用する（スクリープロットで決定）．PCA は線形であるため，scRNA-seq の複雑な非線形構造（異なる細胞型の関係）の可視化には限界がある．[^14][^13]

### UMAP：位相幾何学的非線形次元削減

UMAP（Uniform Manifold Approximation and Projection）は，データが多様体上に分布するという仮定のもと，ファジー単体複体（fuzzy simplicial complex）によって高次元の位相的構造をモデル化し，それを低次元で再現する．[^13][^14]

アルゴリズムの要点：
1. 高次元空間での隣接グラフ構築（\(k\)-最近傍グラフ）
2. ファジー重みの計算（指数関数的減衰）
3. クロスエントロピーを損失関数として低次元埋め込みを最適化

```python
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.pl.umap(adata)
```

主要パラメータ：
- `n_neighbors`：局所 vs 大域構造のトレードオフを制御（小さいと局所的，大きいと大域的）[^14]
- `min_dist`：埋め込み空間での点間距離（小さいほど密集したクラスター）[^14]

> **数学者への注意**：UMAP の埋め込みは確率的アルゴリズムであり，実行ごとに結果が変わる（`random_state` で固定可能）．また，UMAP の距離は位相的に意味を持つが，ユークリッド的な意味での距離は保存しない．クラスター間の見かけ上の近さや遠さを過度に解釈しないこと．[^13]

***

## モジュール 4：クラスタリングとグラフ分割

### Leiden アルゴリズム

scRNA-seq では，まず細胞間の \(k\)-NN グラフを構築し（モジュール 3 で実施済み），次にグラフ分割によるコミュニティ検出を行う．Leiden アルゴリズム（Louvain の改良版）は，グラフのモジュラリティ[^9]

\[ Q = \frac{1}{2m} \sum_{ij} \left[ A_{ij} - \frac{k_i k_j}{2m} \right] \delta(c_i, c_j) \]

を最大化するようにクラスター割り当てを求める（\(A_{ij}\)：隣接行列，\(k_i\)：次数，\(m\)：総辺数）．[^15]

```python
sc.tl.leiden(adata, resolution=0.5)
sc.pl.umap(adata, color=['leiden'])
```

`resolution` パラメータを大きくすると細かいクラスター，小さくすると粗いクラスターが得られる．このパラメータ選択は実質的に「細胞型の解像度」の問題であり，数学的には一意の正解はない．[^16]

***

## モジュール 5：細胞型アノテーション — マーカー遺伝子と仮説検定

### マーカー遺伝子の探索

各クラスターに特異的に高発現する遺伝子（マーカー遺伝子）を統計的に同定する．

```python
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
```

Wilcoxon 順位和検定を用いてクラスター内 vs クラスター外の発現量を比較する．これはノンパラメトリック検定であり，正規性の仮定が不要な点が scRNA-seq の歪んだ分布に適している．[^9]

### 既知マーカー遺伝子によるアノテーション

PBMC の代表的なマーカー遺伝子：[^3][^2]

| 細胞型 | マーカー遺伝子 | 表面マーカー（参考） |
|---|---|---|
| CD4+ T 細胞 | IL7R, CCR7 | CD3+, CD4+ |
| CD8+ T 細胞 | CD8A | CD3+, CD8+ |
| B 細胞 | MS4A1 (CD20), CD79A | CD19+ |
| NK 細胞 | GNLY, NKG7 | CD56+, CD3- |
| 単球（古典的） | CD14, LYZ | CD14+ |
| 単球（非古典的） | FCGR3A, MS4A7 | CD16+ |
| 樹状細胞 | FCER1A, CST3 | — |
| 血小板 | PPBP | — |

```python
marker_genes = ['IL7R', 'CD79A', 'MS4A1', 'CD8A', 'CD14', 'LYZ', 'GNLY', 'FCGR3A', 'FCER1A', 'PPBP']
sc.pl.umap(adata, color=marker_genes)
```

アノテーションを手動で割り当てる：

```python
cell_type_map = {'0': 'CD4 T', '1': 'CD14 Mono', '2': 'B', ...}
adata.obs['cell_type'] = adata.obs['leiden'].map(cell_type_map)
sc.pl.umap(adata, color=['cell_type'])
```

> **自動化ツール**：**ScType** はマーカー遺伝子データベースを用いた全自動のアノテーションツールである．**SingleR** は参照データセットとのスピアマン相関に基づく手法で広く使われる．[^17][^18]

***

## モジュール 6：差次発現解析（DEG） — 負の二項分布モデル

### なぜ負の二項分布か

scRNA-seq の UMI カウントデータは，遺伝子 \(j\) に対してポアソン分布の混合（過分散）として近似され，負の二項分布でよくフィットされる：[^19][^20]

\[ P(X = k) = \binom{k + r - 1}{k} p^k (1-p)^r \]

ここで \(r\) は分散パラメータ（形状パラメータ），\(p\) は成功確率である．平均 \(\mu\) と分散 \(\sigma^2 = \mu + \mu^2/r\) の関係が二次的であることが，RNA-seq データの過分散を表現する上で重要である．[^20]

### 典型的な比較：細胞型間

```python
# 特定の 2 クラスター間の DEG
sc.tl.rank_genes_groups(adata, groupby='cell_type', groups=['CD14 Mono'], 
                         reference='CD4 T', method='wilcoxon')
sc.pl.rank_genes_groups(adata, groups=['CD14 Mono'])
```

**擬似バルク法（pseudobulk）** は同一サンプルの複数細胞を集計してからサンプル間比較を行う手法で，偽陽性を抑制する現在のベストプラクティスである．[^21]

> **統計的注意**：同一被験者由来の多数の細胞を独立サンプルとして扱うと偽陽性が激増する（擬似反復問題）．ランダム効果モデルまたは被験者ごとの集計（pseudobulk）による対処が推奨される．[^21]

***

## モジュール 7：発展 — 軌跡推定と擬似時間

### PAGA（Partition-based Graph Abstraction）

細胞状態の連続的変化（分化過程など）を推定する手法．クラスター間のコネクティビティ（グラフの辺の重み）を計算し，細胞状態空間の位相的構造を可視化する．[^22][^23]

```python
sc.tl.paga(adata, groups='leiden')
sc.pl.paga(adata, plot=True)
```

### 拡散擬似時間（Diffusion Pseudotime, DPT）

拡散過程（マルコフ連鎖）に基づき，参照細胞（根細胞）からの擬似的な「発生時間」を計算する：[^22]

```python
adata.uns['iroot'] = np.flatnonzero(adata.obs['cell_type'] == 'CD34+')
sc.tl.dpt(adata)
sc.pl.umap(adata, color=['dpt_pseudotime'])
```

### RNA velocity（参考）

スプライシング前後の mRNA 比率から細胞の発現状態の「速度ベクトル」を推定する手法が **scVelo** である．これは常微分方程式系のパラメータ推定問題として定式化できる：[^23]

\[ \frac{du}{dt} = \alpha(t) - \beta u, \quad \frac{ds}{dt} = \beta u - \gamma s \]

ここで \(u\) は未スプライシング mRNA，\(s\) はスプライシング済み mRNA，\(\alpha\) は転写率，\(\beta\) はスプライシング率，\(\gamma\) は分解率である．[^23]

***

## ツールと教材のまとめ

| カテゴリ | ツール・リソース | 言語・形式 |
|---|---|---|
| メイン解析 | [Scanpy](https://scanpy.readthedocs.io/) | Python |
| データ形式 | AnnData | Python |
| クラスタリング | Leiden (leidenalg) | Python |
| 次元削減 | umap-learn | Python |
| アンビエント RNA 除去 | SoupX | R |
| ダブレット検出 | scDblFinder | R/Python |
| 自動アノテーション | ScType, SingleR | R/Python |
| 軌跡推定 | scVelo, CellRank | Python |
| 公開データ | [10x PBMC 3k](https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz) | mtx 形式 |
| オンライン教科書 | [sc-best-practices.org](https://www.sc-best-practices.org/) | Web/Jupyter |
| 公式チュートリアル | [Scanpy PBMC 3k](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html) | Jupyter |
| 参考論文 | Luecken & Theis (2019) *Mol Syst Biol* | 査読済み論文 |

***

## 数学者向け概念マッピング

scRNA-seq 解析の各ステップと数理的概念の対応を整理する：

| 解析ステップ | 対応する数理概念 |
|---|---|
| カウント行列 | 非負整数値疎行列 \(\mathbb{Z}_{\geq 0}^{m \times n}\) |
| QC フィルタリング | 外れ値検出（ロバスト統計） |
| 正規化 | スケーリング変換，対数変換 |
| HVG 選択 | 分散の正則化推定，特徴選択 |
| PCA | 特異値分解（SVD），主成分分析 |
| \(k\)-NN グラフ | 計量空間上のグラフ構築 |
| UMAP | 位相的データ解析，多様体学習 |
| Leiden クラスタリング | グラフのモジュラリティ最適化 |
| Wilcoxon 検定 | ノンパラメトリック仮説検定 |
| 負の二項分布モデル | 過分散カウントデータの確率モデル |
| 擬似時間 | マルコフ連鎖，拡散過程 |
| RNA velocity | 常微分方程式系，パラメータ推定 |

***

## 推奨する学習順序と自習リソース

1. **Scanpy 公式チュートリアル**（PBMC 3k）：コードをそのまま実行し，各ステップの出力を確認する[^5][^4]
2. **Galaxy Training Network の Scanpy チュートリアル**：ブラウザ上で GUI ベースに実行できるため，環境構築が不要[^24]
3. **sc-best-practices.org**：各ステップのベストプラクティスと数理的背景の詳細[^6]
4. **Luecken & Theis (2019)**：scRNA-seq 解析のベストプラクティス論文（*Mol Syst Biol*）[^25][^7]
5. **AJACS 講習会動画**（統合 TV, 2025）：日本語による scRNA-seq 解析入門（聴講推奨）[^26]
6. **Dimensionality Reduction in scRNA-Seq（PDF）**：PCA，t-SNE，UMAP の数理的解説[^14]

---

## References

1. [How single cell sequencing data analysis works - 10x Genomics](https://www.10xgenomics.com/blog/how-single-cell-sequencing-data-analysis-works) - Editor's note: The purpose of this piece is to introduce how scRNA-seq data analysis works. Some ste...

2. [Peripheral Blood | Whole blood | Handbook | Miltenyi Biotec | Belgique](https://www.miltenyibiotec.com/BE-en/support/macs-handbook/human-cells-and-organs/human-cell-sources/blood-human.html) - Common lymphoid progenitor cells give rise to T and B lymphocytes and natural killer (NK) cells, whe...

3. [Peripheral Blood Mononuclear Cells - NCBI - NIH](https://www.ncbi.nlm.nih.gov/books/NBK500157/) - PBMCs include lymphocytes (T cells, B cells, and NK cells), monocytes, and dendritic cells. In human...

4. [Preprocessing and clustering 3k PBMCs (legacy workflow) - Scanpy](https://scanpy.readthedocs.io/en/1.10.x/tutorials/basics/clustering-2017.html) - Integrating spatial data with scRNA-seq using scanorama · Experimental · How to preprocess UMI count...

5. [Preprocessing and clustering 3k PBMCs (legacy workflow)](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html) - In May 2017, this started out as a demonstration that Scanpy would allow to reproduce most of Seurat...

6. [Single-Cell Analysis: Best Practice Guidelines - Helmholtz Munich](https://www.helmholtz-munich.de/computational-health-center/news/article/single-cell-analysis-best-practice-guidelines) - Fabian Theis from Helmholtz Munich and the Technical University of Munich (TUM), have developed Cell...

7. [Current best practices in single‐cell RNA‐seq analysis: a tutorial](https://pmc.ncbi.nlm.nih.gov/articles/PMC6582955/) - A detailed case study that integrates the established current best practices is available on our git...

8. [1. Scanpy: Preprocessing and clustering 3k PBMCs - Read the Docs](https://singlecellanalysistutorial.readthedocs.io/en/latest/notebooks/Scanpy_PBMC.html) - Scanpyを用いたクラスタリング解析の基本的なワークフローを紹介します。 Google ColabまたはJupyter notebook上で作業を行います。 内容はSeuratのGuided tut...

9. [[PDF] Single-cell RNA-seq analysis - BaRC](http://barc.wi.mit.edu/education/hot_topics/scRNAseq_2025/SingleCellRNAseq2025.pdf) - Seurat is an R package designed for QC, analysis, and exploration of single cell RNA-seq data. • Dev...

10. [Best Practices For Single-Cell Analysis Across Modalities - Scribd](https://www.scribd.com/document/784728880/Best-practices-for-single-cell-analysis-across-modalities) - Best practices for single-cell analysis across modalities - Free download as PDF File (.pdf), Text F...

11. [SoupX removes ambient RNA contamination from droplet-based ...](https://pmc.ncbi.nlm.nih.gov/articles/PMC7763177/) - We present SoupX, a tool for removing ambient RNA contamination from droplet-based single-cell RNA s...

12. [Best Practices for Analysis of 10x Genomics Single Cell RNA-seq Data](https://www.10xgenomics.com/analysis-guides/best-practices-analysis-10x-single-cell-rnaseq-data) - In this guide, we provide an overview of best practices for analyzing single cell gene expression da...

13. [9. Dimensionality Reduction - Single-cell best practices](https://www.sc-best-practices.org/preprocessing_visualization/dimensionality_reduction.html) - As previously mentioned, scRNA-seq is a high-throughput sequencing technology that produces datasets...

14. [[PDF] Dimensionality Reduction in single-cell RNA-seq PCA, tSNE, UMAP](https://bix.unil.ch/ref/Dimension_Reduction_PCA_TSNE_UMAP.pdf) - Dimensionality reduction of scRNA-Seq data. Gene. Description. Cell 1. Cell 2 ... PCA in scRNA-Seq d...

15. [Analyzing scRNA-seq data by CCP-assisted UMAP and tSNE - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC11642954/) - The aim of this study is to explore CCP's potential as the primary dimensionality reduction method f...

16. [W20: Single-Cell RNA-Seq Analysis with R](https://qcb.ucla.edu/collaboratory/workshops/w20-single-cell-rna-seq/) - This workshop aims to introduce the basic concepts and algorithms for single-cell RNA-seq analysis. ...

17. [Automatic Cell Type Annotation Using Marker Genes for ... - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC9599378/) - The reference scRNA-seq-data-based cell type annotation methods can be divided into several modellin...

18. [Fully-automated and ultra-fast cell-type identification using specific ...](https://www.nature.com/articles/s41467-022-28803-w) - Using six scRNA-seq datasets from various human and mouse tissues, we show how ScType provides unbia...

19. [Evaluating negative binomial models and parameter estimations in ...](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1014014) - Single-cell RNA sequencing (scRNA-seq) measures mRNA molecule counts in individual cells. For most g...

20. [negative bionomial distribution in (single-cell) RNAseq - Chatomics](https://divingintogeneticsandgenomics.com/post/negative-bionomial-distribution-in-single-cell-rnaseq/) - It seems that for each gene, the counts across all cells in scRNAseq data can be modeled with negati...

21. [16. Differential gene expression analysis - Single-cell best practices](https://www.sc-best-practices.org/conditions/differential_gene_expression.html) - Single-cell specific methods were found to be especially prone to wrongly labeling highly expressed ...

22. [Trajectory inference using PAGA - GitHub Pages](https://nbisweden.github.io/workshop-scRNAseq/labs/scanpy/scanpy_07_trajectory.html) - Cellrank can use different metrics to calculate a transition matrix, such as RNA velocity or pseudot...

23. [RNA Velocity Basics — scVelo documentation - Read the Docs](https://scvelo.readthedocs.io/en/stable/VelocityBasics.html) - Here you will learn the basics of RNA velocity analysis. For illustration, it is applied to endocrin...

24. [Clustering 3K PBMCs with Scanpy - Galaxy Training!](https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/scrna-scanpy-pbmc3k/tutorial.html) - In this tutorial, we will investigate clustering of single-cell data from 10x Genomics, including pr...

25. [Current best practices in single-cell RNA-seq analysis: a tutorial](https://pubmed.ncbi.nlm.nih.gov/31217225/) - Here, we detail the steps of a typical single-cell RNA-seq analysis, including pre-processing (quali...

26. [AJACS「シングルセルRNA-seqを知って・学んで・使う」 - YouTube](https://www.youtube.com/watch?v=nathCDFZpEg) - https://doi.org/10.7875/togotv.2025.008 本日の統合TVは、2024年12月23日に開催されたデータ解析講習会：AJACS「シングルセルRNA-seqを知って・学...

